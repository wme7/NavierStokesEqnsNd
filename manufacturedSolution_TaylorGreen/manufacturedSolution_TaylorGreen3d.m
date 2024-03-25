function [q, qx, qy, qz, qxx, qyy, qzz, qt] = manufacturedSolution_TaylorGreen3d(x, y, z, L, time)

global mu

% Refs.:
% [1] Antuono, M. (2020). Tri-periodic fully three-dimensional analytic
%     solutions for the Navier–Stokes equations. Journal of Fluid
%     Mechanics, 890, A23.

% Parameters
mu = 0.025;
rho = 1.0;
nu = mu / rho;

% Model Parameters
A = 4 * sqrt(2) / (3 * sqrt(3));
c56 = 5/6 * pi;
c16 = 1/6 * pi;
p0 = 0.0;
model = 2;
switch model
    case 1, a=c56; b=c16;
    case 2, a=c16; b=c56;
end

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
u = A * (sin(k*xp-a).* cos(k*yp-b).* sin(k*zp) - cos(k*zp-a).* sin(k*xp-b).* sin(k*yp)) * F(time);
v = A * (sin(k*yp-a).* cos(k*zp-b).* sin(k*xp) - cos(k*xp-a).* sin(k*yp-b).* sin(k*zp)) * F(time);
w = A * (sin(k*zp-a).* cos(k*xp-b).* sin(k*yp) - cos(k*yp-a).* sin(k*zp-b).* sin(k*xp)) * F(time);
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

rxx = zeros(size(x));
uxx = zeros(size(x));
uxy = zeros(size(x));
uxz = zeros(size(x));
uyy = zeros(size(x));
uzz = zeros(size(x));
pxx = zeros(size(x));

ryy = zeros(size(x));
vxx = zeros(size(x));
vxy = zeros(size(x));
vyy = zeros(size(x));
vyz = zeros(size(x));
vzz = zeros(size(x));
pyy = zeros(size(y));

rzz = zeros(size(x));
wxx = zeros(size(x));
wyy = zeros(size(x));
wxz = zeros(size(x));
wyz = zeros(size(x));
wzz = zeros(size(x));
pzz = zeros(size(z));

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
qxx = [rxx, uxx, uxy, uxz, uyy, uzz, pxx];
qyy = [ryy, vxx, vxy, vyy, vyz, vzz, pyy];
qzz = [rzz, wxx, wyy, wxz, wyz, wzz, pzz];
qt = [rt, ut, vt, wt, pt];

end % funtion