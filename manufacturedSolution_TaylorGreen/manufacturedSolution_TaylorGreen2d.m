function [q, qx, qy, qt] = manufacturedSolution_TaylorGreen2d(x, y, L, time)
    
% Refs.:
% [1] Antuono, M. (2020). Tri-periodic fully three-dimensional analytic
%     solutions for the Navier–Stokes equations. Journal of Fluid
%     Mechanics, 890, A23.

% Parameters
mu = 0.025;
rho = 1.0;
nu = mu / rho;

% Wave number
k = 2 * pi / L;

% Define exponential function
F = @(t) exp(-2 * nu * k*k * t);

% Exact solution
r = rho * ones(size(x));
u =-cos(k*x).* sin(k*y) * F(time);
v = sin(k*x).* cos(k*y) * F(time);
p =-0.25 * rho * (cos(2*k*x) + cos(2*k*y)) * F(time)^2;

% Compute derivatives
rx = zeros(size(x));
ux = zeros(size(x));
vx = zeros(size(x));
px = zeros(size(x));

ry = zeros(size(y));
uy = zeros(size(y));
vy = zeros(size(y));
py = zeros(size(y));

rt = zeros(size(y));
ut = zeros(size(y));
vt = zeros(size(y));
pt = zeros(size(y));

% Store results in output arrays
q  = [r,  u,  v,  p];
qx = [rx, ux, vx, px];
qy = [ry, uy, vy, py];
qt = [rt, ut, vt, pt];

end % funtion