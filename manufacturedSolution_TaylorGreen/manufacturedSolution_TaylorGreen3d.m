function [q, qx, qy, qz, qt] = manufacturedSolution_TaylorGreen3d(x, y, z, L, time)

% Refs.:
% [1] Antuono, M. (2020). Tri-periodic fully three-dimensional analytic
%     solutions for the Navier–Stokes equations. Journal of Fluid
%     Mechanics, 890, A23.

% Parameters
mu = 0.025;
rho = 1.0;
nu = mu / rho;

% Model Parameters
model = 2;
A = 4 * sqrt(2) / (3 * sqrt(3));
c56 = 5/6 * pi;
c16 = 1/6 * pi;
p0 = 0.0;

% Wave number
k = 2 * pi / L;

% Primed coordinates
R = 0;
psi = acos(R / sqrt(1 + R));
xp = x + psi/k * 1;
yp = y + psi/k * 1;
zp = z + psi/k * 1;

% Define exponential function
F = @(t) exp(-3 * nu * k^2 * t);

% Exact solution
r = rho * ones(size(x));
switch model
    case 1
    u = A * (sin(k*xp-c56).* cos(k*yp-c16).* sin(k*zp) - cos(k*zp-c56).* sin(k*xp-c16).* sin(k*yp)) * F(time);
    v = A * (sin(k*yp-c56).* cos(k*zp-c16).* sin(k*xp) - cos(k*xp-c56).* sin(k*yp-c16).* sin(k*zp)) * F(time);
    w = A * (sin(k*zp-c56).* cos(k*xp-c16).* sin(k*yp) - cos(k*yp-c56).* sin(k*zp-c16).* sin(k*xp)) * F(time);
    case 2
    u = A * (sin(k*xp-c16).* cos(k*yp-c56).* sin(k*zp) - cos(k*zp-c16).* sin(k*xp-c56).* sin(k*yp)) * F(time);
    v = A * (sin(k*yp-c16).* cos(k*zp-c56).* sin(k*xp) - cos(k*xp-c16).* sin(k*yp-c56).* sin(k*zp)) * F(time);
    w = A * (sin(k*zp-c16).* cos(k*xp-c56).* sin(k*yp) - cos(k*yp-c16).* sin(k*zp-c56).* sin(k*xp)) * F(time);
end
p = p0 - rho.*(0.5*(u.*u + v.*v + w.*w) + 0);

% Compute derivatives
rx = zeros(size(x));
ux = zeros(size(x));
vx = zeros(size(x));
wx = zeros(size(x));
px = zeros(size(x));

ry = zeros(size(x));
uy = zeros(size(y));
vy = zeros(size(y));
wy = zeros(size(y));
py = zeros(size(y));

rz = zeros(size(x));
uz = zeros(size(z));
vz = zeros(size(z));
wz = zeros(size(z));
pz = zeros(size(z));

rt = zeros(size(x));
ut = zeros(size(y));
vt = zeros(size(y));
wt = zeros(size(y));
pt = zeros(size(y));

% Store results in output arrays
q = [r, u, v, w, p];
qx = [rx, ux, vx, wx, px];
qy = [ry, uy, vy, wy, py];
qz = [rz, uz, vz, wz, pz];
qt = [rt, ut, vt, wt, pt];

end % funtion