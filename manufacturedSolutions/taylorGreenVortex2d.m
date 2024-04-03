function [q, qt, qx, qy, qxx, qyy, source] = taylorGreenVortex2d(x, y, L, t)

global gamma mu Pr %#ok<GVMIS>

% Refs.:
% [1] Antuono, M. (2020). Tri-periodic fully three-dimensional analytic
%     solutions for the Navier–Stokes equations. Journal of Fluid
%     Mechanics, 890, A23.

% Model parameters
p0 = 1.0;

% Wave number
k = 2 * pi / L;

%% Solution
r = ones(size(x));
u = -exp(-2*k.^2.*mu.*t).*sin(k.*y).*cos(k.*x);
v = exp(-2*k.^2.*mu.*t).*sin(k.*x).*cos(k.*y);
p = p0 - (0.25*cos(2*k.*x) + 0.25*cos(2*k.*y)).*exp(-4*k.^2.*mu.*t);

% Derivatives
r_t = zeros(size(x));
r_x = zeros(size(x));
r_y = zeros(size(x));
r_xx = zeros(size(x));
r_yy = zeros(size(x));
u_t = 2*k.^2.*mu.*exp(-2*k.^2.*mu.*t).*sin(k.*y).*cos(k.*x);
u_x = k.*exp(-2*k.^2.*mu.*t).*sin(k.*x).*sin(k.*y);
u_y = -k.*exp(-2*k.^2.*mu.*t).*cos(k.*x).*cos(k.*y);
u_xx = k.^2.*exp(-2*k.^2.*mu.*t).*sin(k.*y).*cos(k.*x);
u_xy = k.^2.*exp(-2*k.^2.*mu.*t).*sin(k.*x).*cos(k.*y);
u_yy = k.^2.*exp(-2*k.^2.*mu.*t).*sin(k.*y).*cos(k.*x);
v_t = -2*k.^2.*mu.*exp(-2*k.^2.*mu.*t).*sin(k.*x).*cos(k.*y);
v_x = k.*exp(-2*k.^2.*mu.*t).*cos(k.*x).*cos(k.*y);
v_y = -k.*exp(-2*k.^2.*mu.*t).*sin(k.*x).*sin(k.*y);
v_xx = -k.^2.*exp(-2*k.^2.*mu.*t).*sin(k.*x).*cos(k.*y);
v_xy = -k.^2.*exp(-2*k.^2.*mu.*t).*sin(k.*y).*cos(k.*x);
v_yy = -k.^2.*exp(-2*k.^2.*mu.*t).*sin(k.*x).*cos(k.*y);
p_t = 4*k.^2.*mu.*(0.25*cos(2*k.*x) + 0.25*cos(2*k.*y)).*exp(-4*k.^2.*mu.*t);
p_x = 0.5*k.*exp(-4*k.^2.*mu.*t).*sin(2*k.*x);
p_y = 0.5*k.*exp(-4*k.^2.*mu.*t).*sin(2*k.*y);
p_xx = 1.0*k.^2.*exp(-4*k.^2.*mu.*t).*cos(2*k.*x);
p_yy = 1.0*k.^2.*exp(-4*k.^2.*mu.*t).*cos(2*k.*y);

% Store results in output arrays
q  = [r,  u,  v,  p];
qt = [r_t, u_t, v_t, p_t];
qx = [r_x, u_x, v_x, p_x];
qy = [r_y, u_y, v_y, p_y];
qxx = [r_xx, u_xx, u_xy, u_yy, p_xx];
qyy = [r_yy, v_xx, v_xy, v_yy, p_yy];

%% Source Terms
Q_rho = zeros(size(x));
Q_rhou = zeros(size(x));
Q_rhov = zeros(size(x));
Q_rhoE = k.*(-4.0*Pr.*gamma.*k.*mu.*exp(2*k.^2.*mu.*t).*sin(k.*x).^2.*sin(k.*y).^2 + 4.0*Pr.*k.*mu.*exp(2*k.^2.*mu.*t).*sin(k.*x).^2.*sin(k.*y).^2 - 2.0*Pr.*k.*mu.*exp(2*k.^2.*mu.*t).*sin(k.*x).^2 - 2.0*Pr.*k.*mu.*exp(2*k.^2.*mu.*t).*sin(k.*y).^2 + 2.0*Pr.*k.*mu.*exp(2*k.^2.*mu.*t) + 1.0*Pr.*sin(k.*x).^3.*sin(k.*y) - 1.0*Pr.*sin(k.*x).*sin(k.*y).^3 + 2.0*gamma.*k.*mu.*exp(2*k.^2.*mu.*t).*sin(k.*x).^2 + 2.0*gamma.*k.*mu.*exp(2*k.^2.*mu.*t).*sin(k.*y).^2 - 2.0*gamma.*k.*mu.*exp(2*k.^2.*mu.*t)).*exp(-6*k.^2.*mu.*t)./(Pr.*(gamma - 1));

source = [Q_rho, Q_rhou, Q_rhov, Q_rhoE];
end % funtion