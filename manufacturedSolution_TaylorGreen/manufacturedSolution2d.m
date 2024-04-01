function [q, qt, qx, qy, qxx, qyy, source] = manufacturedSolution2d(x, y, L, t)

% Refs.:
% [1] Malaya, N., Estacio-Hiroms, K. C., Stogner, R. H., Schulz, K. W.,
%     Bauman, P. T., & Carey, G. F. (2013). MASA: a library for
%     verification using manufactured and analytical solutions. Engineering
%     with Computers, 29, 487-496.

global gamma mu Pr %#ok<*GVMIS>

% Parameters
r_0 = 1.0; r_T = 0.1; r_X = 0.1; r_Y = 0.2;
u_0 = 0.0; u_T = 0.1; u_X = 0.1; u_Y = 0.2;
v_0 = 0.0; v_T = 0.1; v_X = 0.1; v_Y = 0.2;
p_0 = 1.0; p_T = 0.1; p_X = 0.1; p_Y = 0.2;

artPiInvL = 2*pi/L; arxPiInvL = 2*pi/L; aryPiInvL = 2*pi/L;
autPiInvL = 2*pi/L; auxPiInvL = 2*pi/L; auyPiInvL = 2*pi/L;
avtPiInvL = 2*pi/L; avxPiInvL = 2*pi/L; avyPiInvL = 2*pi/L;
aptPiInvL = 2*pi/L; apxPiInvL = 2*pi/L; apyPiInvL = 2*pi/L;

%% Solution
r = r_0 + r_T.*sin(artPiInvL.*t) + r_X.*sin(arxPiInvL.*x) + r_Y.*cos(aryPiInvL.*y);
u = u_0 + u_T.*cos(autPiInvL.*t) + u_X.*sin(auxPiInvL.*x) + u_Y.*cos(auyPiInvL.*y);
v = v_0 + v_T.*sin(avtPiInvL.*t) + v_X.*cos(avxPiInvL.*x) + v_Y.*sin(avyPiInvL.*y);
p = p_0 + p_T.*cos(aptPiInvL.*t) + p_X.*cos(apxPiInvL.*x) + p_Y.*sin(apyPiInvL.*y);

% Derivatives
r_t = artPiInvL.*r_T.*cos(artPiInvL.*t);
r_x = arxPiInvL.*r_X.*cos(arxPiInvL.*x);
r_y = -aryPiInvL.*r_Y.*sin(aryPiInvL.*y);
r_xx = -arxPiInvL.^2.*r_X.*sin(arxPiInvL.*x);
r_yy = -aryPiInvL.^2.*r_Y.*cos(aryPiInvL.*y);
u_t = -autPiInvL.*u_T.*sin(autPiInvL.*t);
u_x = auxPiInvL.*u_X.*cos(auxPiInvL.*x);
u_y = -auyPiInvL.*u_Y.*sin(auyPiInvL.*y);
u_xx = -auxPiInvL.^2.*u_X.*sin(auxPiInvL.*x);
u_xy = zeros(size(x));
u_yy = -auyPiInvL.^2.*u_Y.*cos(auyPiInvL.*y);
v_t = avtPiInvL.*v_T.*cos(avtPiInvL.*t);
v_x = -avxPiInvL.*v_X.*sin(avxPiInvL.*x);
v_y = avyPiInvL.*v_Y.*cos(avyPiInvL.*y);
v_xx = -avxPiInvL.^2.*v_X.*cos(avxPiInvL.*x);
v_xy = zeros(size(x));
v_yy = -avyPiInvL.^2.*v_Y.*sin(avyPiInvL.*y);
p_t = -aptPiInvL.*p_T.*sin(aptPiInvL.*t);
p_x = -apxPiInvL.*p_X.*sin(apxPiInvL.*x);
p_y = apyPiInvL.*p_Y.*cos(apyPiInvL.*y);
p_xx = -apxPiInvL.^2.*p_X.*cos(apxPiInvL.*x);
p_yy = -apyPiInvL.^2.*p_Y.*sin(apyPiInvL.*y);

% Store results in output arrays
q  = [r,  u,  v,  p];
qt = [r_t, u_t, v_t, p_t];
qx = [r_x, u_x, v_x, p_x];
qy = [r_y, u_y, v_y, p_y];
qxx = [r_xx, u_xx, u_xy, u_yy, p_xx];
qyy = [r_yy, v_xx, v_xy, v_yy, p_yy];

%% Source Terms & auxiliary variables
RHO_2 = r;
U_2 = u;
V_2 = v;
P_2 = p;

% Source Terms
Q_rho = RHO_2.*auxPiInvL.*u_X.*cos(auxPiInvL.*x) + RHO_2.*avyPiInvL.*v_Y.*cos(avyPiInvL.*y) + U_2.*arxPiInvL.*r_X.*cos(arxPiInvL.*x) - V_2.*aryPiInvL.*r_Y.*sin(aryPiInvL.*y) + artPiInvL.*r_T.*cos(artPiInvL.*t);
Q_rhou = 2*RHO_2.*U_2.*auxPiInvL.*u_X.*cos(auxPiInvL.*x) + RHO_2.*U_2.*avyPiInvL.*v_Y.*cos(avyPiInvL.*y) - RHO_2.*V_2.*auyPiInvL.*u_Y.*sin(auyPiInvL.*y) - RHO_2.*autPiInvL.*u_T.*sin(autPiInvL.*t) + U_2.^2.*arxPiInvL.*r_X.*cos(arxPiInvL.*x) - U_2.*V_2.*aryPiInvL.*r_Y.*sin(aryPiInvL.*y) + U_2.*artPiInvL.*r_T.*cos(artPiInvL.*t) - apxPiInvL.*p_X.*sin(apxPiInvL.*x) + 1.33333333333333*auxPiInvL.^2.*mu.*u_X.*sin(auxPiInvL.*x) + auyPiInvL.^2.*mu.*u_Y.*cos(auyPiInvL.*y);
Q_rhov = -RHO_2.*U_2.*avxPiInvL.*v_X.*sin(avxPiInvL.*x) + RHO_2.*V_2.*auxPiInvL.*u_X.*cos(auxPiInvL.*x) + 2*RHO_2.*V_2.*avyPiInvL.*v_Y.*cos(avyPiInvL.*y) + RHO_2.*avtPiInvL.*v_T.*cos(avtPiInvL.*t) + U_2.*V_2.*arxPiInvL.*r_X.*cos(arxPiInvL.*x) - V_2.^2.*aryPiInvL.*r_Y.*sin(aryPiInvL.*y) + V_2.*artPiInvL.*r_T.*cos(artPiInvL.*t) + apyPiInvL.*p_Y.*cos(apyPiInvL.*y) + avxPiInvL.^2.*mu.*v_X.*cos(avxPiInvL.*x) + 1.33333333333333*avyPiInvL.^2.*mu.*v_Y.*sin(avyPiInvL.*y);
Q_rhoE = P_2.*auxPiInvL.*u_X.*cos(auxPiInvL.*x) + P_2.*avyPiInvL.*v_Y.*cos(avyPiInvL.*y) + RHO_2.*U_2.*(U_2.*auxPiInvL.*u_X.*cos(auxPiInvL.*x) - V_2.*avxPiInvL.*v_X.*sin(avxPiInvL.*x) + (-P_2.*arxPiInvL.*r_X.*cos(arxPiInvL.*x) - RHO_2.*apxPiInvL.*p_X.*sin(apxPiInvL.*x))./(RHO_2.^2.*(gamma - 1))) + RHO_2.*V_2.*(-U_2.*auyPiInvL.*u_Y.*sin(auyPiInvL.*y) + V_2.*avyPiInvL.*v_Y.*cos(avyPiInvL.*y) + (P_2.*aryPiInvL.*r_Y.*sin(aryPiInvL.*y) + RHO_2.*apyPiInvL.*p_Y.*cos(apyPiInvL.*y))./(RHO_2.^2.*(gamma - 1))) + RHO_2.*auxPiInvL.*u_X.*(P_2./(RHO_2.*(gamma - 1)) + U_2.^2/2 + V_2.^2/2).*cos(auxPiInvL.*x) + RHO_2.*avyPiInvL.*v_Y.*(P_2./(RHO_2.*(gamma - 1)) + U_2.^2/2 + V_2.^2/2).*cos(avyPiInvL.*y) + RHO_2.*(-U_2.*autPiInvL.*u_T.*sin(autPiInvL.*t) + V_2.*avtPiInvL.*v_T.*cos(avtPiInvL.*t) + (-P_2.*artPiInvL.*r_T.*cos(artPiInvL.*t) - RHO_2.*aptPiInvL.*p_T.*sin(aptPiInvL.*t))./(RHO_2.^2.*(gamma - 1))) - U_2.*apxPiInvL.*p_X.*sin(apxPiInvL.*x) + U_2.*arxPiInvL.*r_X.*(P_2./(RHO_2.*(gamma - 1)) + U_2.^2/2 + V_2.^2/2).*cos(arxPiInvL.*x) + 1.33333333333333*U_2.*auxPiInvL.^2.*mu.*u_X.*sin(auxPiInvL.*x) + U_2.*auyPiInvL.^2.*mu.*u_Y.*cos(auyPiInvL.*y) + V_2.*apyPiInvL.*p_Y.*cos(apyPiInvL.*y) - V_2.*aryPiInvL.*r_Y.*(P_2./(RHO_2.*(gamma - 1)) + U_2.^2/2 + V_2.^2/2).*sin(aryPiInvL.*y) + V_2.*avxPiInvL.^2.*mu.*v_X.*cos(avxPiInvL.*x) + 1.33333333333333*V_2.*avyPiInvL.^2.*mu.*v_Y.*sin(avyPiInvL.*y) + artPiInvL.*r_T.*(P_2./(RHO_2.*(gamma - 1)) + U_2.^2/2 + V_2.^2/2).*cos(artPiInvL.*t) - auxPiInvL.*u_X.*(2*auxPiInvL.*mu.*u_X.*cos(auxPiInvL.*x) - 0.666666666666667*mu.*(auxPiInvL.*u_X.*cos(auxPiInvL.*x) + avyPiInvL.*v_Y.*cos(avyPiInvL.*y))).*cos(auxPiInvL.*x) + auyPiInvL.*mu.*u_Y.*(-auyPiInvL.*u_Y.*sin(auyPiInvL.*y) - avxPiInvL.*v_X.*sin(avxPiInvL.*x)).*sin(auyPiInvL.*y) + avxPiInvL.*mu.*v_X.*(-auyPiInvL.*u_Y.*sin(auyPiInvL.*y) - avxPiInvL.*v_X.*sin(avxPiInvL.*x)).*sin(avxPiInvL.*x) - avyPiInvL.*v_Y.*(2*avyPiInvL.*mu.*v_Y.*cos(avyPiInvL.*y) - 0.666666666666667*mu.*(auxPiInvL.*u_X.*cos(auxPiInvL.*x) + avyPiInvL.*v_Y.*cos(avyPiInvL.*y))).*cos(avyPiInvL.*y) - gamma.*mu.*(-RHO_2.^2.*apxPiInvL.^2.*p_X.*cos(apxPiInvL.*x) - RHO_2.*(-P_2.*arxPiInvL.^2.*r_X.*sin(arxPiInvL.*x) - 2*apxPiInvL.*arxPiInvL.*p_X.*r_X.*sin(apxPiInvL.*x).*cos(arxPiInvL.*x)) + arxPiInvL.^2.*r_X.^2.*(2*p_0 + 2*p_T.*cos(aptPiInvL.*t) + 2*p_X.*cos(apxPiInvL.*x) + 2*p_Y.*sin(apyPiInvL.*y)).*cos(arxPiInvL.*x).^2)./(Pr.*RHO_2.^3.*(gamma - 1)) - gamma.*mu.*(-RHO_2.^2.*apyPiInvL.^2.*p_Y.*sin(apyPiInvL.*y) - RHO_2.*(-P_2.*aryPiInvL.^2.*r_Y.*cos(aryPiInvL.*y) - 2*apyPiInvL.*aryPiInvL.*p_Y.*r_Y.*sin(aryPiInvL.*y).*cos(apyPiInvL.*y)) + aryPiInvL.^2.*r_Y.^2.*(2*p_0 + 2*p_T.*cos(aptPiInvL.*t) + 2*p_X.*cos(apxPiInvL.*x) + 2*p_Y.*sin(apyPiInvL.*y)).*sin(aryPiInvL.*y).^2)./(Pr.*RHO_2.^3.*(gamma - 1));

source = [Q_rho, Q_rhou, Q_rhov, Q_rhoE];
end % funtion