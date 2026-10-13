% SBSP Pointing Simulation - Week 2
% Attitude control: B-dot detumbling + PD nadir pointing
% Same orbit as Week 1 (600 km SSO)
%
% Phase 1 - the satellite comes off the launcher tumbling. B-dot uses the
%           magnetorquers to brake it against Earth's magnetic field.
% Phase 2 - once slow enough, reaction wheels + PD controller lock the
%           z-axis (where the power beam comes out) onto nadir.
%           Magnetorquers switch to dumping wheel momentum in the background.
%
% Quaternion convention: q = [q1 q2 q3 q4], q4 = scalar part.
% A = q_to_dcm(q) maps ECI vectors into the body frame: v_body = A*v_eci
% (Markley & Crassidis, Fundamentals of Spacecraft Attitude Det. & Control)

clear; clc; close all;

%% constants
mu      = 398600.4418;   % km^3/s^2
R_E     = 6378.137;      % km
J2      = 1.08263e-3;
omega_E = 7.2921150e-5;  % rad/s

%% orbit (same as Week 1)
h_alt = 600;                  % km
a     = R_E + h_alt;
inc   = deg2rad(97.81);       % SSO inclination

RAAN_0 = deg2rad(90);
u_0    = 0;

n_orb = sqrt(mu / a^3);       % mean motion [rad/s]
T_orb = 2*pi / n_orb;         % period [s]

RAAN_dot = -(3/2)*n_orb*J2*(R_E/a)^2*cos(inc);
u_dot    = n_orb*(1 + (3/4)*J2*(R_E/a)^2*(5*cos(inc)^2 - 1));

%% simulation time
n_orbits = 8;
dt = 0.5;                     % s  (RK4, plenty small for these rates)
t  = 0:dt:(n_orbits*T_orb);
N  = length(t);

fprintf('--- Week 2: attitude control ---\n')
fprintf('orbits   : %d (%.1f min)\n', n_orbits, t(end)/60)
fprintf('timestep : %.2f s (%d steps)\n\n', dt, N)

%% orbit propagation (ECI, km and km/s)
RAAN = RAAN_0 + RAAN_dot*t;
u    = u_0    + u_dot*t;

Rx = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];

X_eci = zeros(3,N);
V_eci = zeros(3,N);
for k = 1:N
    cO = cos(RAAN(k)); sO = sin(RAAN(k));
    Rz = [cO -sO 0; sO cO 0; 0 0 1];
    R_oi = Rz*Rx;

    X_eci(:,k) = R_oi * (a*[cos(u(k)); sin(u(k)); 0]);
    V_eci(:,k) = R_oi * (a*u_dot*[-sin(u(k)); cos(u(k)); 0]);
end

%% spacecraft inertia
% ~200 kg bus with the solar arrays along y, so Iy is the biggest
J    = diag([45 85 50]);     % kg*m^2
Jinv = inv(J);

%% Earth magnetic field - tilted dipole
% mu0/(4*pi) * Earth's dipole moment (8.0e22 A*m^2)
B0 = 7.94e15;                % T*m^3
% geomagnetic north pole ~80.7 N, 72.7 W (IGRF-13). The dipole moment
% itself points the opposite way, toward the southern hemisphere.
pole_lat = deg2rad(80.7);
pole_lon = deg2rad(-72.7);
m_hat = -[cos(pole_lat)*cos(pole_lon);
          cos(pole_lat)*sin(pole_lon);
          sin(pole_lat)];

%% actuators
m_max   = 10;      % A*m^2   magnetorquer limit, per axis
tau_max = 0.05;    % N*m     reaction wheel torque limit, per axis
h_max   = 1.0;     % N*m*s   wheel momentum capacity (only checked, not enforced)

%% controller gains

% B-dot: m = -k*dB/dt. Gain is high on purpose - the rods sit at their limit
% until the rates are down to ~0.1 deg/s, i.e. it behaves like bang-bang
% B-dot, which is the fastest you can detumble with a given rod size.
k_bdot = 1e8;      % A*m^2 / (T/s)

% PD: tau = -Kp*dq - Kd*w_err
% the quaternion vector part is ~theta/2 for small angles, hence the 2 in Kp
% so that wn and zeta below are the real closed-loop values on each axis
wn   = 0.03;       % rad/s  (~3.5 min period)
zeta = 0.9;
Kp = 2*wn^2*J;
Kd = 2*zeta*wn*J;

w_slew = deg2rad(0.25);   % max slew rate for big maneuvers [rad/s]

% wheel momentum dumping with the magnetorquers (cross-product law)
k_dump = 1e-3;     % 1/s

% handover from B-dot to PD
w_handover = deg2rad(0.5);    % rad/s

%% initial conditions
w = deg2rad([3.5; -2.8; 4.1]);           % tumble after separation [rad/s]

ax0 = [1; 1; 1]/sqrt(3);                 % some arbitrary starting attitude
q   = [ax0*sin(deg2rad(30)/2); cos(deg2rad(30)/2)];

h_rw = zeros(3,1);                       % wheels start at rest

% magnetic field at t = 0, so the first B-dot finite difference is sane
B_prev = q_to_dcm(q) * dipole_field(X_eci(:,1), 0, B0, m_hat, omega_E);

%% storage
w_hist    = zeros(3,N);
h_hist    = zeros(3,N);
tau_mtq_h = zeros(3,N);
tau_rw_h  = zeros(3,N);
err_hist  = zeros(1,N);   % 3-axis attitude error [deg]
mode_hist = zeros(1,N);   % 0 = B-dot, 1 = PD

%% main loop
fprintf('simulating...\n')
mode = 0;
k_ho = [];

for k = 1:N
    r = X_eci(:,k);
    v = V_eci(:,k);
    A = q_to_dcm(q);

    % --- reference attitude: LVLH / nadir pointing ---
    % z -> nadir (beam axis), y -> minus orbit normal, x -> roughly velocity
    h_orb = cross(r, v);
    z_ref = -r/norm(r);
    y_ref = -h_orb/norm(h_orb);
    x_ref = cross(y_ref, z_ref);
    A_ref = [x_ref'; y_ref'; z_ref'];

    % the LVLH frame turns once per orbit, about the orbit normal
    w_ref_eci = h_orb / dot(r, r);

    % --- attitude error, in the body frame ---
    dq = dcm_to_quat(A * A_ref');
    if dq(4) < 0
        dq = -dq;            % q and -q are the same attitude, take the short way
    end
    w_err = w - A*w_ref_eci;

    err_hist(k)  = rad2deg(2*acos(min(1, dq(4))));
    w_hist(:,k)  = w;
    h_hist(:,k)  = h_rw;
    mode_hist(k) = mode;

    if k == N
        break
    end

    % --- environment ---
    B_body = A * dipole_field(r, t(k), B0, m_hat, omega_E);

    % gravity gradient, the main disturbance torque in LEO
    r_b    = A*(r/norm(r));
    tau_gg = 3*mu/norm(r)^3 * cross(r_b, J*r_b);

    % --- control ---
    if mode == 0
        % B-dot detumbling, magnetorquers only
        B_dot = (B_body - B_prev)/dt;
        m_cmd = sat(-k_bdot*B_dot, m_max);

        tau_mtq = cross(m_cmd, B_body);
        tau_rw  = zeros(3,1);

        if norm(w) < w_handover && t(k) > 60
            mode = 1;
            k_ho = k;
            fprintf('  B-dot -> PD at t = %.1f min (orbit %.2f)\n', t(k)/60, t(k)/T_orb)
        end
    else
        % PD on the reaction wheels. Same as tau = -Kp*dq - Kd*w_err, just
        % written as a rate command (Kd\Kp = wn/zeta) so the slew rate can be
        % capped. Without the cap a big slew after handover spins the body
        % up to ~1 deg/s and dumps ~2 N*m*s into the wheels.
        w_cmd = -(wn/zeta)*dq(1:3);
        if norm(w_cmd) > w_slew
            w_cmd = w_cmd*w_slew/norm(w_cmd);
        end
        % + gyroscopic feed-forward so it doesn't fight the w x (Jw + h) coupling
        tau_rw = -Kd*(w_err - w_cmd) + cross(w, J*w + h_rw);
        tau_rw = sat(tau_rw, tau_max);

        % dump wheel momentum: m = k (h x B)/|B|^2  gives  m x B = -k*h_perp
        % this torque goes into the body, the PD loop holds the attitude,
        % so the wheels end up absorbing it and slowing down.
        m_cmd   = sat(k_dump*cross(h_rw, B_body)/dot(B_body, B_body), m_max);
        tau_mtq = cross(m_cmd, B_body);
    end
    B_prev = B_body;

    tau_mtq_h(:,k) = tau_mtq;
    tau_rw_h(:,k)  = tau_rw;

    % --- propagate attitude with RK4 ---
    % state x = [w; q; h_rw], torques held constant over the step
    tau_ext = tau_mtq + tau_gg;
    x  = [w; q; h_rw];
    f1 = att_rhs(x,           tau_ext, tau_rw, J, Jinv);
    f2 = att_rhs(x + dt/2*f1, tau_ext, tau_rw, J, Jinv);
    f3 = att_rhs(x + dt/2*f2, tau_ext, tau_rw, J, Jinv);
    f4 = att_rhs(x + dt*f3,   tau_ext, tau_rw, J, Jinv);
    x  = x + dt/6*(f1 + 2*f2 + 2*f3 + f4);

    w    = x(1:3);
    q    = x(4:7)/norm(x(4:7));
    h_rw = x(8:10);
end

%% results
fprintf('\n--- results ---\n')
if isempty(k_ho)
    fprintf('never left B-dot mode - try more orbits or a bigger k_bdot\n')
else
    k_fine = find(err_hist(k_ho:end) < 1, 1) + k_ho - 1;
    if ~isempty(k_fine)
        fprintf('error < 1 deg after    : %.1f min from handover\n', (t(k_fine)-t(k_ho))/60)
    end
    k_ss = find(t > t(k_ho) + T_orb, 1);    % skip the first orbit after handover
    if ~isempty(k_ss)
        fprintf('steady-state error     : %.3f deg max, %.3f deg rms\n', ...
            max(err_hist(k_ss:end)), sqrt(mean(err_hist(k_ss:end).^2)))
    end
end
fprintf('final error            : %.3f deg\n', err_hist(end))
fprintf('max wheel momentum     : %.3f N*m*s (capacity %.1f)\n', max(abs(h_hist(:))), h_max)
fprintf('max wheel torque       : %.1f mN*m\n\n', 1e3*max(abs(tau_rw_h(:))))

%% plots
t_min = t/60;
c_bg   = [0.05 0.05 0.12];
c_x    = [1.00 0.35 0.30];
c_y    = [0.30 0.90 0.40];
c_z    = [0.30 0.55 1.00];
c_grid = [0.25 0.25 0.35];
c_txt  = [0.70 0.70 0.80];
c_orng = [1.00 0.55 0.10];
c_err  = [0.90 0.65 0.15];

figure('Name','Week 2 - Attitude Control','Color',c_bg,'Position',[50 50 1500 950]);

% 1) body rates
ax1 = subplot(4,1,1); style_ax(ax1, c_bg, c_grid, c_txt);
p1 = plot(t_min, rad2deg(w_hist(1,:)),'Color',c_x,'LineWidth',1.2);
p2 = plot(t_min, rad2deg(w_hist(2,:)),'Color',c_y,'LineWidth',1.2);
p3 = plot(t_min, rad2deg(w_hist(3,:)),'Color',c_z,'LineWidth',1.2);
ylabel('\omega [deg/s]','Color','w');
title('Detumbling + Nadir Pointing','Color','w','FontSize',14,'FontWeight','bold');
mark_handover(t_min, k_ho, c_orng, true);
legend([p1 p2 p3], {'\omega_x','\omega_y','\omega_z'}, 'TextColor','w', ...
    'Color',[0.08 0.08 0.15],'EdgeColor',[0.4 0.4 0.5],'Location','northeast');

% 2) pointing error, log scale so the sub-degree part is actually visible
ax2 = subplot(4,1,2); style_ax(ax2, c_bg, c_grid, c_txt);
plot(t_min, max(err_hist, 1e-5),'Color',c_err,'LineWidth',1.3);
set(ax2,'YScale','log');
ylim(ax2,[1e-5 300]);
yline(1,'--','Color',c_orng,'LineWidth',1,'HandleVisibility','off', ...
    'Label','1 deg','LabelColor',c_orng,'LabelHorizontalAlignment','left');
ylabel('Attitude error [deg]','Color','w');
mark_handover(t_min, k_ho, c_orng, false);

% 3) actuator torque magnitudes. log scale again - the handover slew is
%    ~50 mN*m while B-dot and steady pointing are 100-1000x smaller
tq_rw  = 1e3*vecnorm(tau_rw_h);
tq_mtq = 1e3*vecnorm(tau_mtq_h);
tq_rw(tq_rw == 0)   = NaN;      % wheels are off during B-dot
tq_mtq(tq_mtq == 0) = NaN;

ax3 = subplot(4,1,3); style_ax(ax3, c_bg, c_grid, c_txt);
p1 = plot(t_min, tq_mtq,'Color',[0.85 0.85 0.85],'LineWidth',1.0);
p2 = plot(t_min, tq_rw, 'Color',c_z,'LineWidth',1.0);
set(ax3,'YScale','log');
ylim(ax3,[1e-6 1e2]);
ylabel('|Torque| [mN\cdotm]','Color','w');
mark_handover(t_min, k_ho, c_orng, false);
legend([p1 p2], {'magnetorquers','reaction wheels'}, 'TextColor','w', ...
    'Color',[0.08 0.08 0.15],'EdgeColor',[0.4 0.4 0.5],'Location','northeast');

% 4) wheel momentum - shows the dumping working
ax4 = subplot(4,1,4); style_ax(ax4, c_bg, c_grid, c_txt);
p1 = plot(t_min, h_hist(1,:),'Color',c_x,'LineWidth',1.2);
p2 = plot(t_min, h_hist(2,:),'Color',c_y,'LineWidth',1.2);
p3 = plot(t_min, h_hist(3,:),'Color',c_z,'LineWidth',1.2);
ylabel('Wheel momentum [N\cdotm\cdots]','Color','w');
xlabel('Time [min]','Color','w');
mark_handover(t_min, k_ho, c_orng, false);
legend([p1 p2 p3], {'h_x','h_y','h_z'}, 'TextColor','w', ...
    'Color',[0.08 0.08 0.15],'EdgeColor',[0.4 0.4 0.5],'Location','northeast');

linkaxes([ax1 ax2 ax3 ax4],'x');
xlim(ax1,[0 t_min(end)]);


%% ---------------- local functions ----------------

function dx = att_rhs(x, tau_ext, tau_rw, J, Jinv)
% rigid body + reaction wheels
%   J*w_dot = tau_ext + tau_rw - w x (J*w + h)
%   h_dot   = -tau_rw   (wheel spins up opposite to the torque it puts on the body)
%   q_dot   = 0.5 * [q4*w - w x qv ; -w.qv]
w  = x(1:3);
qv = x(4:6);
q4 = x(7);
h  = x(8:10);

w_dot = Jinv*(tau_ext + tau_rw - cross(w, J*w + h));
q_dot = 0.5*[q4*w - cross(w, qv); -dot(w, qv)];
h_dot = -tau_rw;

dx = [w_dot; q_dot; h_dot];
end

function B_eci = dipole_field(r_eci_km, t, B0, m_hat, omega_E)
% tilted dipole, field in tesla, returned in ECI
% (Greenwich assumed aligned with the ECI x-axis at t = 0)
th = omega_E*t;
R_ie = [cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1];   % ECI -> ECEF

r = R_ie*r_eci_km*1e3;               % m, ECEF
rn = norm(r);
rh = r/rn;
B_ecef = B0/rn^3 * (3*dot(m_hat, rh)*rh - m_hat);

B_eci = R_ie'*B_ecef;
end

function y = sat(x, lim)
y = min(max(x, -lim), lim);
end

function A = q_to_dcm(q)
% q = [q1 q2 q3 q4], q4 scalar  ->  A such that v_body = A*v_eci
q1=q(1); q2=q(2); q3=q(3); q4=q(4);
A = [q4^2+q1^2-q2^2-q3^2,  2*(q1*q2+q3*q4),      2*(q1*q3-q2*q4);
     2*(q1*q2-q3*q4),      q4^2-q1^2+q2^2-q3^2,  2*(q2*q3+q1*q4);
     2*(q1*q3+q2*q4),      2*(q2*q3-q1*q4),      q4^2-q1^2-q2^2+q3^2];
end

function q = dcm_to_quat(A)
% inverse of q_to_dcm (Shepperd's method - use the biggest of the four
% candidates as the divisor so nothing blows up near 180 deg)
tr = trace(A);
[~, i] = max([A(1,1) A(2,2) A(3,3) tr]);
switch i
    case 4
        s = 2*sqrt(1 + tr);
        q = [(A(2,3)-A(3,2))/s; (A(3,1)-A(1,3))/s; (A(1,2)-A(2,1))/s; s/4];
    case 1
        s = 2*sqrt(1 + A(1,1) - A(2,2) - A(3,3));
        q = [s/4; (A(1,2)+A(2,1))/s; (A(1,3)+A(3,1))/s; (A(2,3)-A(3,2))/s];
    case 2
        s = 2*sqrt(1 + A(2,2) - A(1,1) - A(3,3));
        q = [(A(1,2)+A(2,1))/s; s/4; (A(2,3)+A(3,2))/s; (A(3,1)-A(1,3))/s];
    case 3
        s = 2*sqrt(1 + A(3,3) - A(1,1) - A(2,2));
        q = [(A(1,3)+A(3,1))/s; (A(2,3)+A(3,2))/s; s/4; (A(1,2)-A(2,1))/s];
end
q = q/norm(q);
end

function style_ax(ax, c_bg, c_grid, c_txt)
set(ax,'Color',c_bg,'XColor',c_txt,'YColor',c_txt,'GridColor',c_grid,'GridAlpha',0.4);
hold(ax,'on'); grid(ax,'on');
end

function mark_handover(t_min, k_ho, c, with_label)
if isempty(k_ho), return, end
if with_label
    xline(t_min(k_ho),'--','Color',c,'LineWidth',1.2,'HandleVisibility','off', ...
        'Label','PD takeover','LabelColor',c,'LabelVerticalAlignment','top');
else
    xline(t_min(k_ho),'--','Color',c,'LineWidth',1.2,'HandleVisibility','off');
end
end