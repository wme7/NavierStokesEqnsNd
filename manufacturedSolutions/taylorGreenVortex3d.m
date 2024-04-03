function [q, qt, qx, qy, qz, qxx, qyy, qzz, source] = taylorGreenVortex3d(x, y, z, L, t)

global gamma mu Pr %#ok<GVMIS>

% Refs.:
% [1] Antuono, M. (2020). Tri-periodic fully three-dimensional analytic
%     solutions for the Navier–Stokes equations. Journal of Fluid
%     Mechanics, 890, A23.

% Model parameters
p0 = 1.0; [A, B, C] = deal(0.5);

% Wave number
k = 2 * pi / L;

%% Solution
r = ones(size(x));
u = (A.*sin(k.*z) + C.*cos(k.*y)).*exp(-k.^2.*mu.*t);
v = (A.*cos(k.*z) + B.*sin(k.*x)).*exp(-k.^2.*mu.*t);
w = (B.*cos(k.*x) + C.*sin(k.*y)).*exp(-k.^2.*mu.*t);
p = p0 - (A.*B.*sin(k.*x).*cos(k.*z) + A.*C.*sin(k.*z).*cos(k.*y) + B.*C.*sin(k.*y).*cos(k.*x)).*exp(-2*k.^2.*mu.*t);

% Derivatives
r_t = zeros(size(x));
r_x = zeros(size(x));
r_y = zeros(size(x));
r_z = zeros(size(x));
r_xx = zeros(size(x));
r_yy = zeros(size(x));
r_zz = zeros(size(x));
u_t = -k.^2.*mu.*(A.*sin(k.*z) + C.*cos(k.*y)).*exp(-k.^2.*mu.*t);
u_x = zeros(size(x));
u_y = -C.*k.*exp(-k.^2.*mu.*t).*sin(k.*y);
u_z = A.*k.*exp(-k.^2.*mu.*t).*cos(k.*z);
u_xx = zeros(size(x));
u_xy = zeros(size(x));
u_xz = zeros(size(x));
u_yy = -C.*k.^2.*exp(-k.^2.*mu.*t).*cos(k.*y);
u_zz = -A.*k.^2.*exp(-k.^2.*mu.*t).*sin(k.*z);
v_t = -k.^2.*mu.*(A.*cos(k.*z) + B.*sin(k.*x)).*exp(-k.^2.*mu.*t);
v_x = B.*k.*exp(-k.^2.*mu.*t).*cos(k.*x);
v_y = zeros(size(x));
v_z = -A.*k.*exp(-k.^2.*mu.*t).*sin(k.*z);
v_xx = -B.*k.^2.*exp(-k.^2.*mu.*t).*sin(k.*x);
v_xy = zeros(size(x));
v_yy = zeros(size(x));
v_yz = zeros(size(x));
v_zz = -A.*k.^2.*exp(-k.^2.*mu.*t).*cos(k.*z);
w_t = -k.^2.*mu.*(B.*cos(k.*x) + C.*sin(k.*y)).*exp(-k.^2.*mu.*t);
w_x = -B.*k.*exp(-k.^2.*mu.*t).*sin(k.*x);
w_y = C.*k.*exp(-k.^2.*mu.*t).*cos(k.*y);
w_z = zeros(size(x));
w_xx = -B.*k.^2.*exp(-k.^2.*mu.*t).*cos(k.*x);
w_xz = zeros(size(x));
w_yy = -C.*k.^2.*exp(-k.^2.*mu.*t).*sin(k.*y);
w_yz = zeros(size(x));
w_zz = zeros(size(x));
p_t = 2*k.^2.*mu.*(A.*B.*sin(k.*x).*cos(k.*z) + A.*C.*sin(k.*z).*cos(k.*y) + B.*C.*sin(k.*y).*cos(k.*x)).*exp(-2*k.^2.*mu.*t);
p_x = -(A.*B.*k.*cos(k.*x).*cos(k.*z) - B.*C.*k.*sin(k.*x).*sin(k.*y)).*exp(-2*k.^2.*mu.*t);
p_y = -(-A.*C.*k.*sin(k.*y).*sin(k.*z) + B.*C.*k.*cos(k.*x).*cos(k.*y)).*exp(-2*k.^2.*mu.*t);
p_z = -(-A.*B.*k.*sin(k.*x).*sin(k.*z) + A.*C.*k.*cos(k.*y).*cos(k.*z)).*exp(-2*k.^2.*mu.*t);
p_xx = -(-A.*B.*k.^2.*sin(k.*x).*cos(k.*z) - B.*C.*k.^2.*sin(k.*y).*cos(k.*x)).*exp(-2*k.^2.*mu.*t);
p_yy = -(-A.*C.*k.^2.*sin(k.*z).*cos(k.*y) - B.*C.*k.^2.*sin(k.*y).*cos(k.*x)).*exp(-2*k.^2.*mu.*t);
p_zz = -(-A.*B.*k.^2.*sin(k.*x).*cos(k.*z) - A.*C.*k.^2.*sin(k.*z).*cos(k.*y)).*exp(-2*k.^2.*mu.*t);

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
Q_rhoE = k.*(A.^2.*B.*Pr.*sin(k.*(x - 2*z))/4 - A.^2.*B.*Pr.*sin(k.*(x + 2*z))/4 + A.^2.*C.*Pr.*cos(k.*(y - 2*z))/4 - A.^2.*C.*Pr.*cos(k.*(y + 2*z))/4 - A.^2.*Pr.*gamma.*k.*mu.*exp(k.^2.*mu.*t) + A.^2.*Pr.*k.*mu.*exp(k.^2.*mu.*t) + A.*B.^2.*Pr.*cos(k.*(2*x - z))/4 - A.*B.^2.*Pr.*cos(k.*(2*x + z))/4 + 3*A.*B.*C.*Pr.*sin(k.*(-x + y + z))/4 + 3*A.*B.*C.*Pr.*sin(k.*(x - y + z))/4 + 3*A.*B.*C.*Pr.*sin(k.*(x + y - z))/4 - 3*A.*B.*C.*Pr.*sin(k.*(x + y + z))/4 - 3*A.*B.*C.*Pr.*cos(k.*(-x + y + z))/4 - 3*A.*B.*C.*Pr.*cos(k.*(x - y + z))/4 - 3*A.*B.*C.*Pr.*cos(k.*(x + y - z))/4 - 3*A.*B.*C.*Pr.*cos(k.*(x + y + z))/4 + A.*B.*Pr.*gamma.*k.*mu.*exp(k.^2.*mu.*t).*sin(k.*(x - z)) + A.*B.*Pr.*gamma.*k.*mu.*exp(k.^2.*mu.*t).*sin(k.*(x + z)) - A.*B.*gamma.*k.*mu.*exp(k.^2.*mu.*t).*sin(k.*(x - z)) - A.*B.*gamma.*k.*mu.*exp(k.^2.*mu.*t).*sin(k.*(x + z)) - A.*C.^2.*Pr.*sin(k.*(2*y - z))/4 - A.*C.^2.*Pr.*sin(k.*(2*y + z))/4 - A.*C.*Pr.*gamma.*k.*mu.*exp(k.^2.*mu.*t).*sin(k.*(y - z)) + A.*C.*Pr.*gamma.*k.*mu.*exp(k.^2.*mu.*t).*sin(k.*(y + z)) + A.*C.*gamma.*k.*mu.*exp(k.^2.*mu.*t).*sin(k.*(y - z)) - A.*C.*gamma.*k.*mu.*exp(k.^2.*mu.*t).*sin(k.*(y + z)) - B.^2.*C.*Pr.*sin(k.*(2*x - y))/4 - B.^2.*C.*Pr.*sin(k.*(2*x + y))/4 - B.^2.*Pr.*gamma.*k.*mu.*exp(k.^2.*mu.*t) + B.^2.*Pr.*k.*mu.*exp(k.^2.*mu.*t) + B.*C.^2.*Pr.*cos(k.*(x - 2*y))/4 - B.*C.^2.*Pr.*cos(k.*(x + 2*y))/4 - B.*C.*Pr.*gamma.*k.*mu.*exp(k.^2.*mu.*t).*sin(k.*(x - y)) + B.*C.*Pr.*gamma.*k.*mu.*exp(k.^2.*mu.*t).*sin(k.*(x + y)) + B.*C.*gamma.*k.*mu.*exp(k.^2.*mu.*t).*sin(k.*(x - y)) - B.*C.*gamma.*k.*mu.*exp(k.^2.*mu.*t).*sin(k.*(x + y)) - C.^2.*Pr.*gamma.*k.*mu.*exp(k.^2.*mu.*t) + C.^2.*Pr.*k.*mu.*exp(k.^2.*mu.*t)).*exp(-3*k.^2.*mu.*t)./(Pr.*(gamma - 1));

source = [Q_rho, Q_rhou, Q_rhov, Q_rhow, Q_rhoE];
end % funtion