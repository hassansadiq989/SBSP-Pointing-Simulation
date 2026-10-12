% SBSP Pointing Simulation - Week 1
% Orbital propagator + sun vector model for a sun-synchronous orbit
% Rectenna at Kirchberg Park, Luxembourg (49.6233N, 6.1557E)
%
% Run with F5. Left = spinning Earth + satellite + beam, right = elevation angle.
% Earth rotates under a fixed orbit, so the rectenna stays over Luxembourg.
% Sun position lets us detect eclipses (no power when in Earth's shadow).
% TODO: tune satellite size / beam angle once the attitude controller is in

clear; clc; close all;

%% constants
% standard values (Vallado, 2013)
mu      = 398600.4418;   % km^3/s^2
R_E     = 6378.137;      % km
J2      = 1.08263e-3;    % oblateness - needed for SSO
omega_E = 7.2921150e-5;  % rad/s, Earth spin rate

%% orbit parameters
% 600 km SSO - 97.81 deg inclination gives the right nodal regression
h = 600;
a = R_E + h;
i = deg2rad(97.81);

RAAN_0 = deg2rad(90);   % orbit plane orientation at t=0
u_0    = 0;             % start at ascending node

%% derived quantities
n = sqrt(mu / a^3);   % mean motion [rad/s]
T = 2*pi / n;         % period [s], ~97 min

% J2 nodal regression - makes it sun-synchronous (need +0.9856 deg/day)
RAAN_dot = -(3/2) * n * J2 * (R_E/a)^2 * cos(i);
u_dot    = n * (1 + (3/4)*J2*(R_E/a)^2*(5*cos(i)^2 - 1));

fprintf('--- orbit check ---\n')
fprintf('altitude   : %d km\n', h)
fprintf('period     : %.1f min\n', T/60)
fprintf('RAAN drift : %+.4f deg/day  (want +0.9856)\n', rad2deg(RAAN_dot)*86400)
fprintf('\n')

%% time vector
n_orbits = 5;
N  = n_orbits * 600;
t  = linspace(0, n_orbits*T, N);
dt = t(2) - t(1);

%% propagate orbit in ECI
% satellite pos = circle in orbital plane, rotated by inclination then RAAN
RAAN = RAAN_0 + RAAN_dot * t;
u    = u_0    + u_dot    * t;

% inclination rotation is constant - compute once
Rx = [1,      0,       0   ;
      0,  cos(i),  -sin(i) ;
      0,  sin(i),   cos(i)];

X_eci = zeros(3, N);
for k = 1:N
    r_orb = a * [cos(u(k)); sin(u(k)); 0];
    cO = cos(RAAN(k));  sO = sin(RAAN(k));
    Rz = [cO -sO 0; sO cO 0; 0 0 1];
    X_eci(:,k) = Rz * Rx * r_orb;
end

%% rectenna position (Kirchberg Park)
% keep it in ECEF (fixed to the ground). The Earth group gets rotated in
% the animation, so the marker stays on Luxembourg automatically.
lat_r = deg2rad(49.6233);
lon_r = deg2rad(6.1557);
r_rect_ecef = R_E * [cos(lat_r)*cos(lon_r);
                     cos(lat_r)*sin(lon_r);
                     sin(lat_r)           ];

% current ECI position of the rectenna at each step (for geometry + beam)
r_rect_eci = zeros(3, N);
for k = 1:N
    th = omega_E * t(k);
    Rz_E = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
    r_rect_eci(:,k) = Rz_E * r_rect_ecef;
end

%% elevation angle + access windows
elev_min = 5;
elev   = zeros(1, N);
access = false(1, N);
for k = 1:N
    dr = X_eci(:,k) - r_rect_eci(:,k);
    up = r_rect_eci(:,k) / norm(r_rect_eci(:,k));
    elev(k)   = rad2deg(asin(dot(dr/norm(dr), up)));
    access(k) = elev(k) > elev_min;
end

n_passes     = sum(diff([false, access]) == 1);
total_access = sum(access) * dt;
fprintf('--- Kirchberg Park visibility ---\n')
fprintf('passes (5 orbits) : %d\n', n_passes)
fprintf('total access time : %.1f min\n', total_access/60)
if n_passes > 0
    fprintf('avg pass duration : %.1f min\n', (total_access/60)/n_passes)
end
fprintf('\n')

%% sun position in ECI
% approximate algorithm from Vallado (2013) ch.5 - good to ~0.01 deg
% we need a calendar date to know where the sun is. pick an epoch:
epoch_JD = 2460676.5;   % Jan 1 2025, 00:00 UT1 (a recent round number)

sun_eci = zeros(3, N);
in_sun  = true(1, N);   % true = sunlit, false = eclipse

for k = 1:N
    % Julian centuries since J2000.0
    JD_now = epoch_JD + t(k)/86400;
    T_ut1  = (JD_now - 2451545.0) / 36525;

    % mean longitude + mean anomaly of the sun [deg]
    lam_M = mod(280.460 + 36000.771*T_ut1, 360);
    M_sun = mod(357.5277 + 35999.050*T_ut1, 360);

    % ecliptic longitude (add equation of center)
    lam_ec = lam_M + 1.9147*sind(M_sun) + 0.0200*sind(2*M_sun);

    % obliquity of the ecliptic
    eps = 23.4393 - 0.0130*T_ut1;

    % sun distance [AU] - not critical, we just need the direction
    r_AU = 1.00014 - 0.01671*cosd(M_sun) - 0.00014*cosd(2*M_sun);

    % unit vector from Earth to Sun, in ECI (equatorial frame)
    sun_eci(:,k) = [cosd(lam_ec);
                    cosd(eps)*sind(lam_ec);
                    sind(eps)*sind(lam_ec)];

    % --- eclipse check (cylindrical shadow model) ---
    % project satellite position onto the Earth-Sun direction.
    % if the projection is negative (sat is on the dark side) AND the
    % perpendicular distance from the Earth-Sun axis is < R_E, we're in shadow.
    s_hat = sun_eci(:,k);
    proj  = dot(X_eci(:,k), s_hat);          % along sun direction
    perp  = norm(X_eci(:,k) - proj*s_hat);   % sideways from the axis
    if proj < 0 && perp < R_E
        in_sun(k) = false;
    end
end

% power can only be beamed when we're (a) in sunlight and (b) over the rectenna
can_beam = access & in_sun;

n_eclipse  = sum(~in_sun);
pct_sunlit = 100*sum(in_sun)/N;
beam_time  = sum(can_beam)*dt;

fprintf('--- sun & eclipse ---\n')
fprintf('epoch          : JD %.1f\n', epoch_JD)
fprintf('sunlit fraction: %.1f%%\n', pct_sunlit)
fprintf('eclipse time   : %.1f min (of %.1f total)\n', n_eclipse*dt/60, t(end)/60)
fprintf('beam-able time : %.1f min (access AND sunlit)\n', beam_time/60)
fprintf('\n')

%% colors
c_bg    = [0.05 0.05 0.12];
% (Earth color comes from topo texture now, no need for manual RGB)
c_orbit = [0.35 0.75 1.00];
c_green = [0.15 1.00 0.45];   % access arcs
c_red   = [1.00 0.30 0.20];   % beam searching / rectenna
c_blue  = [0.20 0.65 1.00];   % beam locked on rectenna
c_body  = [0.85 0.85 0.88];   % satellite body
c_panel = [0.10 0.15 0.45];   % solar panels
c_gold  = [0.80 0.70 0.35];   % dish
c_grid  = [0.25 0.25 0.35];
c_txt   = [0.70 0.70 0.80];
c_orng  = [1.00 0.55 0.10];
c_sun   = [1.00 0.85 0.20];   % sun arrow
c_eclp  = [0.50 0.20 0.55];   % eclipse shading on the plot

sat_size = 130;   % km per model unit - exaggerated so the sat is visible

%% figure
fig = figure('Name','SBSP orbit','Color',c_bg,'Position',[60 50 1500 740]);

% ---------- left panel: 3D scene ----------
ax1 = subplot(1,2,1);
set(ax1,'Color',c_bg); hold on; axis equal; axis vis3d; grid on;
ax1.GridColor = c_grid; ax1.GridAlpha = 0.4;
ax1.XColor = c_txt; ax1.YColor = c_txt; ax1.ZColor = c_txt;

% everything glued to the Earth goes under this transform so it can spin
t_earth = hgtransform('Parent',ax1);

% textured Earth using MATLAB's built-in topography data
% topo.mat ships with every MATLAB install - 180x360 elevation grid
load('topo','topo','topomap1');

% build an RGB texture from the elevation data + colormap
% topomap1 is a 64-color map: blues for ocean, greens/browns for land
cmap = topomap1;
topo_norm = topo - min(topo(:));
topo_norm = topo_norm / max(topo_norm(:));                  % 0..1
idx = round(topo_norm * (size(cmap,1)-1)) + 1;             % colormap indices
earth_rgb = zeros(size(topo,1), size(topo,2), 3);
for ch = 1:3
    tmp = cmap(idx, ch);
    earth_rgb(:,:,ch) = reshape(tmp, size(topo));
end

% build our own sphere mesh from explicit lat/lon so the texture aligns
% topo: rows = lat -90..+89 (south→north), cols = lon 0..359
% we want the mesh to match that exactly, then no shifting is needed
lat_grid = linspace(-pi/2, pi/2, size(topo,1));   % south→north like topo
lon_grid = linspace(0, 2*pi, size(topo,2));        % 0→360 like topo
[LON, LAT] = meshgrid(lon_grid, lat_grid);
xs = cos(LAT) .* cos(LON);
ys = cos(LAT) .* sin(LON);
zs = sin(LAT);

surf(R_E*xs, R_E*ys, R_E*zs, 'Parent',t_earth, ...
    'FaceColor','texturemap','CData',earth_rgb, ...
    'EdgeColor','none','FaceAlpha',1.0, ...
    'AmbientStrength',0.45,'DiffuseStrength',0.65,'SpecularStrength',0.15);

% rectenna marker (sits just above the surface, rides with the Earth)
rr = r_rect_ecef*1.012;
plot3(rr(1),rr(2),rr(3),'p','Parent',t_earth, ...
    'Color',c_red,'MarkerFaceColor',c_red,'MarkerSize',16);

% equator ring (ECI, faint reference)
th_eq = linspace(0,2*pi,200);
plot3(R_E*cos(th_eq), R_E*sin(th_eq), zeros(1,200), ...
    'Color',[c_orbit 0.25],'LineWidth',0.8);

% orbit track + green access arcs (fixed in ECI)
h_orbit = plot3(X_eci(1,:),X_eci(2,:),X_eci(3,:),'Color',c_orbit,'LineWidth',1.2);
X_acc = X_eci; X_acc(:,~access) = NaN;
h_access = plot3(X_acc(1,:),X_acc(2,:),X_acc(3,:),'Color',c_green,'LineWidth',3.5);

% --- build the satellite under its own transform (geometry once) ---
t_sat = hgtransform('Parent',ax1);
box_patch(t_sat,[0 0 0],   [0.7 0.7 1.0], c_body, 1);     % body
box_patch(t_sat,[0 0.85 0],[0.08 0.9 0.08], c_body,1);    % boom +y
box_patch(t_sat,[0 -0.85 0],[0.08 0.9 0.08],c_body,1);    % boom -y
box_patch(t_sat,[0 2.0 0], [1.1 2.2 0.05], c_panel,1);    % panel +y
box_patch(t_sat,[0 -2.0 0],[1.1 2.2 0.05], c_panel,1);    % panel -y

% transmit dish on the +z (nadir) face, pointing at Earth
[thd,rd] = meshgrid(linspace(0,2*pi,24), linspace(0,0.45,6));
surf(rd.*cos(thd), rd.*sin(thd), 0.55+0.6*rd.^2, 'Parent',t_sat, ...
    'FaceColor',c_gold,'EdgeColor','none');

% four thruster nozzles on the -z (anti-nadir) face
for sx = [-0.22 0.22]
    for sy = [-0.22 0.22]
        [xn,yn,zn] = cone_surf([sx sy -0.55],[0 0 -1],28,0.28,12);
        surf(xn,yn,zn,'Parent',t_sat,'FaceColor',[0.3 0.3 0.33],'EdgeColor','none');
    end
end

% beam cone (surface we recolor each frame) + a bright link line on lock
[xb,yb,zb] = cone_surf(X_eci(:,1),-X_eci(:,1)/norm(X_eci(:,1)),8,h,20);
h_beam = surf(xb,yb,zb,'FaceColor',c_red,'EdgeColor','none','FaceAlpha',0.22);
h_link = plot3([NaN NaN],[NaN NaN],[NaN NaN],'-','Color',c_blue,'LineWidth',2.5);

% sun direction arrow - big yellow arrow pointing where the sun is
% we'll update it each frame, but draw the initial one here
sun0 = sun_eci(:,1);
sun_len = R_E * 2.5;   % long enough to be clearly visible past the Earth
h_sun_arrow = quiver3(0,0,0, sun0(1)*sun_len, sun0(2)*sun_len, sun0(3)*sun_len, 0, ...
    'Color',c_sun,'LineWidth',2.5,'MaxHeadSize',0.4);

% small sun disc at the tip of the arrow (makes it look nicer)
h_sun_dot = plot3(sun0(1)*sun_len, sun0(2)*sun_len, sun0(3)*sun_len, 'o', ...
    'MarkerFaceColor',c_sun,'MarkerEdgeColor','none','MarkerSize',14);

% proxies so the legend can show satellite + beam + sun
hp_sat  = plot3(NaN,NaN,NaN,'s','MarkerFaceColor',c_panel,'MarkerEdgeColor','none','MarkerSize',9);
hp_beam = plot3(NaN,NaN,NaN,'-','Color',c_blue,'LineWidth',3);
hp_sun  = plot3(NaN,NaN,NaN,'-','Color',c_sun,'LineWidth',3);

xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('SBSP — SSO Orbit (600 km)','Color','w','FontSize',14,'FontWeight','bold');
legend([h_orbit h_access hp_sat hp_beam hp_sun], ...
    {'Orbit','Access (>5 deg)','Satellite','Power beam','Sun direction'}, ...
    'TextColor','w','Color',[0.08 0.08 0.15],'EdgeColor',[0.4 0.4 0.5], ...
    'Location','northeast','FontSize',8);
view(35,22);
lighting gouraud;
camlight('headlight');
material([0.6 0.7 0.15 5]);   % less specular = no washed-out glare on the Earth
rotate3d on;                   % let the user drag to rotate the view mid-animation

% ---------- right panel: elevation angle ----------
ax2 = subplot(1,2,2);
set(ax2,'Color',c_bg); hold on; grid on;
ax2.GridColor = c_grid; ax2.GridAlpha = 0.4;
ax2.XColor = c_txt; ax2.YColor = c_txt;
t_min = t/60;

for k = 1:N-1
    % green tint for access windows, purple tint for eclipse
    if access(k)
        patch([t_min(k) t_min(k+1) t_min(k+1) t_min(k)], [-30 -30 95 95], ...
            c_green,'FaceAlpha',0.10,'EdgeColor','none');
    end
    if ~in_sun(k)
        patch([t_min(k) t_min(k+1) t_min(k+1) t_min(k)], [-30 -30 95 95], ...
            c_eclp,'FaceAlpha',0.15,'EdgeColor','none');
    end
end
plot(t_min, elev, 'Color', c_orbit, 'LineWidth', 1.2);
elev_acc = elev; elev_acc(~access) = NaN;
plot(t_min, elev_acc, 'Color', c_green, 'LineWidth', 2.5);
yline(elev_min,'--','Color',c_orng,'LineWidth',1.2, ...
    'Label',sprintf('  min elevation %d deg',elev_min), ...
    'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left','LabelColor',c_orng);
yline(0,'-','Color',[0.5 0.5 0.6],'LineWidth',0.6);
h_cursor = xline(0,'-','Color',[c_gold 0.8],'LineWidth',1.5);

xlabel('Time [min]','Color','w','FontSize',11);
ylabel('Elevation [deg]','Color','w','FontSize',11);
title('Kirchberg Park — visibility & eclipse','Color','w','FontSize',14,'FontWeight','bold');
ylim([-30 95]);

%% animation
fprintf('animating...\n')
for loop = 1:2
    for k = 1:4:N
        if ~ishandle(fig), break; end
        pos = X_eci(:,k);

        % satellite orientation: z points down at Earth, x along velocity
        if k < N, vel = X_eci(:,k+1)-pos; else, vel = pos-X_eci(:,k-1); end
        zb_ = -pos/norm(pos);
        xb_ = vel - (vel'*zb_)*zb_;  xb_ = xb_/norm(xb_);
        yb_ = cross(zb_,xb_);
        M = eye(4);
        M(1:3,1:3) = [xb_ yb_ zb_]*sat_size;
        M(1:3,4)   = pos;
        t_sat.Matrix = M;

        % spin the Earth group to the current sidereal angle
        t_earth.Matrix = makehgtform('zrotate', omega_E*t(k));

        % update sun arrow (sun barely moves over 5 orbits, but still)
        sv = sun_eci(:,k);
        set(h_sun_arrow,'UData',sv(1)*sun_len,'VData',sv(2)*sun_len,'WData',sv(3)*sun_len);
        set(h_sun_dot,'XData',sv(1)*sun_len,'YData',sv(2)*sun_len,'ZData',sv(3)*sun_len);

        % beam cone from sat down to the surface along nadir
        alt = norm(pos) - R_E;
        [xb2,yb2,zb2] = cone_surf(pos, zb_, 8, alt, 20);
        if ~in_sun(k)
            % eclipse - no solar power, beam is off (dim it out)
            set(h_beam,'XData',xb2,'YData',yb2,'ZData',zb2, ...
                'FaceColor',c_eclp,'FaceAlpha',0.10);
            set(h_link,'XData',[NaN NaN],'YData',[NaN NaN],'ZData',[NaN NaN]);
        elseif access(k)
            % sunlit + over rectenna = actively beaming power
            set(h_beam,'XData',xb2,'YData',yb2,'ZData',zb2, ...
                'FaceColor',c_blue,'FaceAlpha',0.22);
            set(h_link,'XData',[pos(1) r_rect_eci(1,k)], ...
                       'YData',[pos(2) r_rect_eci(2,k)], ...
                       'ZData',[pos(3) r_rect_eci(3,k)]);
        else
            % sunlit but not over rectenna - beam searching
            set(h_beam,'XData',xb2,'YData',yb2,'ZData',zb2, ...
                'FaceColor',c_red,'FaceAlpha',0.22);
            set(h_link,'XData',[NaN NaN],'YData',[NaN NaN],'ZData',[NaN NaN]);
        end

        set(h_cursor,'Value',t_min(k));
        drawnow;
        pause(0.01);
    end
end
fprintf('done - next up: attitude control\n')


%% ---- helper functions ----

function box_patch(parent, c, d, col, alpha)
% draw a cuboid centered at c with dimensions d, as a child of parent
dx=d(1)/2; dy=d(2)/2; dz=d(3)/2;
V = [-dx -dy -dz;  dx -dy -dz;  dx  dy -dz; -dx  dy -dz;
     -dx -dy  dz;  dx -dy  dz;  dx  dy  dz; -dx  dy  dz];
V = V + c(:)';
F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
patch('Parent',parent,'Vertices',V,'Faces',F, ...
      'FaceColor',col,'EdgeColor','none','FaceAlpha',alpha);
end

function [X,Y,Z] = cone_surf(apex, axisdir, halfdeg, len, nseg)
% cone mesh: apex at 'apex', opening along 'axisdir', given half-angle & length
ang = linspace(0,2*pi,nseg);
zc  = [0; 1];                       % 0 at apex, 1 at base
rr  = zc * (len*tand(halfdeg));
Xl  = rr*cos(ang);
Yl  = rr*sin(ang);
Zl  = zc*len .* ones(1,nseg);
R   = rot_a_to_b([0;0;1], axisdir(:)/norm(axisdir));
P   = R*[Xl(:)'; Yl(:)'; Zl(:)'] + apex(:);
X = reshape(P(1,:),2,nseg);
Y = reshape(P(2,:),2,nseg);
Z = reshape(P(3,:),2,nseg);
end

function R = rot_a_to_b(a, b)
% rotation matrix that turns unit vector a into unit vector b
v = cross(a,b); s = norm(v); c = dot(a,b);
if s < 1e-9
    if c > 0, R = eye(3); else, R = diag([1 -1 -1]); end
    return
end
vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
R = eye(3) + vx + vx^2*((1-c)/s^2);
end
