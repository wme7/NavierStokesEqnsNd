function [q, qt, qx, qy, qz, qxx, qyy, qzz, source] = taylorGreenVortex3d_AntuonoModel(x, y, z, L, t)

global gamma mu Pr %#ok<RESWD>

% Refs.:
% [1] Antuono, M. (2020). Tri-periodic fully three-dimensional analytic
%     solutions for the Navier–Stokes equations. Journal of Fluid
%     Mechanics, 890, A23.

% Model parameters
A = 4 * sqrt(2) / (3 * sqrt(3));
c56 = 5/6 * pi;
c16 = 1/6 * pi;
model = 2;
switch model
    case 1, a=c56; b=c16;
    case 2, a=c16; b=c56;
end

% Background pressure
p0 = 0.;

% Wave number
k = 2 * pi / L;

% Primed coordinates
R = 0;
psi = acos(R / sqrt(1 + R));
xp = x + psi/k * 1;
yp = y + psi/k * 1;
zp = z + psi/k * 1;

% Define exponential function
F = @(t) exp(-3 * mu * k^2 * t);

% Exact solution
r = ones(size(x));
u = A * (sin(k*xp-a).* cos(k*yp-b).* sin(k*zp) - cos(k*zp-a).* sin(k*xp-b).* sin(k*yp)) * F(t);
v = A * (sin(k*yp-a).* cos(k*zp-b).* sin(k*xp) - cos(k*xp-a).* sin(k*yp-b).* sin(k*zp)) * F(t);
w = A * (sin(k*zp-a).* cos(k*xp-b).* sin(k*yp) - cos(k*yp-a).* sin(k*zp-b).* sin(k*xp)) * F(t);
p = p0 - r .* (u.*u + v.*v + w.*w) / 2;

% Compute derivatives
r_t = zeros(size(x));
u_t = zeros(size(y));
v_t = zeros(size(y));
w_t = zeros(size(y));
p_t = zeros(size(y));

r_x = zeros(size(x));
u_x = zeros(size(x));
v_x = zeros(size(x));
w_x = zeros(size(x));
p_x = zeros(size(x));

r_y = zeros(size(x));
u_y = zeros(size(y));
v_y = zeros(size(y));
w_y = zeros(size(y));
p_y = zeros(size(y));

r_z = zeros(size(x));
u_z = zeros(size(z));
v_z = zeros(size(z));
w_z = zeros(size(z));
p_z = zeros(size(z));

r_xx = zeros(size(x));
u_xx = zeros(size(x));
u_xy = zeros(size(x));
u_xz = zeros(size(x));
u_yy = zeros(size(x));
u_zz = zeros(size(x));
p_xx = zeros(size(x));

r_yy = zeros(size(x));
v_xx = zeros(size(x));
v_xy = zeros(size(x));
v_yy = zeros(size(x));
v_yz = zeros(size(x));
v_zz = zeros(size(x));
p_yy = zeros(size(y));

r_zz = zeros(size(x));
w_xx = zeros(size(x));
w_yy = zeros(size(x));
w_xz = zeros(size(x));
w_yz = zeros(size(x));
w_zz = zeros(size(x));
p_zz = zeros(size(z));

% Store results in output arrays
q = [r, u, v, w, p];
qt = [r_t, u_t, v_t, w_t, p_t];
qx = [r_x, u_x, v_x, w_x, p_x];
qy = [r_y, u_y, v_y, w_y, p_y];
qz = [r_z, u_z, v_z, w_z, p_z];
qxx = [r_xx, u_xx, u_xy, u_xz, u_yy, u_zz, p_xx];
qyy = [r_yy, v_xx, v_xy, v_yz, v_yy, v_zz, p_yy];
qzz = [r_zz, w_xx, w_xz, w_yz, w_yy, w_zz, p_zz];

%% Source Terms
Q_rho = zeros(size(x));
Q_rhou = zeros(size(x));
Q_rhov = zeros(size(x));
Q_rhow = zeros(size(x));
Q_rhoE = zeros(size(x));

source = [Q_rho, Q_rhou, Q_rhov, Q_rhow, Q_rhoE];

end % funtion