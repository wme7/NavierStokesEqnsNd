function [q, qt, qx, qy, qz, qxx, qyy, qzz, source] = manufacturedSolution3d(x, y, z, L, t)

% Refs.:
% [1] Malaya, N., Estacio-Hiroms, K. C., Stogner, R. H., Schulz, K. W.,
%     Bauman, P. T., & Carey, G. F. (2013). MASA: a library for
%     verification using manufactured and analytical solutions. Engineering
%     with Computers, 29, 487-496.

global gamma mu Pr %#ok<*GVMIS>

% Parameters
r_0 = 1.0; r_T = 0.1; r_X = 0.2; r_Y = 0.3; r_Z = 0.4;
u_0 = 0.0; u_T = 0.1; u_X = 0.2; u_Y = 0.3; u_Z = 0.4;
v_0 = 0.0; v_T = 0.1; v_X = 0.2; v_Y = 0.3; v_Z = 0.4;
w_0 = 0.0; w_T = 0.1; w_X = 0.2; w_Y = 0.3; w_Z = 0.4;
p_0 = 1.0; p_T = 0.1; p_X = 0.2; p_Y = 0.3; p_Z = 0.4;

artPiInvL = 2*pi/L; arxPiInvL = 2*pi/L; aryPiInvL = 2*pi/L; arzPiInvL = 2*pi/L;
autPiInvL = 2*pi/L; auxPiInvL = 2*pi/L; auyPiInvL = 2*pi/L; auzPiInvL = 2*pi/L;
avtPiInvL = 2*pi/L; avxPiInvL = 2*pi/L; avyPiInvL = 2*pi/L; avzPiInvL = 2*pi/L;
awtPiInvL = 2*pi/L; awxPiInvL = 2*pi/L; awyPiInvL = 2*pi/L; awzPiInvL = 2*pi/L;
aptPiInvL = 2*pi/L; apxPiInvL = 2*pi/L; apyPiInvL = 2*pi/L; apzPiInvL = 2*pi/L;

%% Solution
r = r_0 + r_T.*sin(artPiInvL.*t) + r_X.*sin(arxPiInvL.*x) + r_Y.*cos(aryPiInvL.*y) + r_Z.*sin(arzPiInvL.*z);
u = u_0 + u_T.*cos(autPiInvL.*t) + u_X.*sin(auxPiInvL.*x) + u_Y.*cos(auyPiInvL.*y) + u_Z.*cos(auzPiInvL.*z);
v = v_0 + v_T.*sin(avtPiInvL.*t) + v_X.*cos(avxPiInvL.*x) + v_Y.*sin(avyPiInvL.*y) + v_Z.*sin(avzPiInvL.*z);
w = w_0 + w_T.*cos(awtPiInvL.*t) + w_X.*sin(awxPiInvL.*x) + w_Y.*sin(awyPiInvL.*y) + w_Z.*cos(awzPiInvL.*z);
p = p_0 + p_T.*cos(aptPiInvL.*t) + p_X.*cos(apxPiInvL.*x) + p_Y.*sin(apyPiInvL.*y) + p_Z.*cos(apzPiInvL.*z);

% Derivatives
r_t = artPiInvL.*r_T.*cos(artPiInvL.*t);
r_x = arxPiInvL.*r_X.*cos(arxPiInvL.*x);
r_y = -aryPiInvL.*r_Y.*sin(aryPiInvL.*y);
r_z = arzPiInvL.*r_Z.*cos(arzPiInvL.*z);
r_xx = -arxPiInvL.^2.*r_X.*sin(arxPiInvL.*x);
r_yy = -aryPiInvL.^2.*r_Y.*cos(aryPiInvL.*y);
r_zz = -arzPiInvL.^2.*r_Z.*sin(arzPiInvL.*z);
u_t = -autPiInvL.*u_T.*sin(autPiInvL.*t);
u_x = auxPiInvL.*u_X.*cos(auxPiInvL.*x);
u_y = -auyPiInvL.*u_Y.*sin(auyPiInvL.*y);
u_z = -auzPiInvL.*u_Z.*sin(auzPiInvL.*z);
u_xx = -auxPiInvL.^2.*u_X.*sin(auxPiInvL.*x);
u_xy = zeros(size(x));
u_xz = zeros(size(x));
u_yy = -auyPiInvL.^2.*u_Y.*cos(auyPiInvL.*y);
u_zz = -auzPiInvL.^2.*u_Z.*cos(auzPiInvL.*z);
v_t = avtPiInvL.*v_T.*cos(avtPiInvL.*t);
v_x = -avxPiInvL.*v_X.*sin(avxPiInvL.*x);
v_y = avyPiInvL.*v_Y.*cos(avyPiInvL.*y);
v_z = avzPiInvL.*v_Z.*cos(avzPiInvL.*z);
v_xx = -avxPiInvL.^2.*v_X.*cos(avxPiInvL.*x);
v_xy = zeros(size(x));
v_yy = -avyPiInvL.^2.*v_Y.*sin(avyPiInvL.*y);
v_yz = zeros(size(x));
v_zz = -avzPiInvL.^2.*v_Z.*sin(avzPiInvL.*z);
w_t = -awtPiInvL.*w_T.*sin(awtPiInvL.*t);
w_x = awxPiInvL.*w_X.*cos(awxPiInvL.*x);
w_y = awyPiInvL.*w_Y.*cos(awyPiInvL.*y);
w_z = -awzPiInvL.*w_Z.*sin(awzPiInvL.*z);
w_xx = -awxPiInvL.^2.*w_X.*sin(awxPiInvL.*x);
w_xz = zeros(size(x));
w_yy = -awyPiInvL.^2.*w_Y.*sin(awyPiInvL.*y);
w_yz = zeros(size(x));
w_zz = -awzPiInvL.^2.*w_Z.*cos(awzPiInvL.*z);
p_t = -aptPiInvL.*p_T.*sin(aptPiInvL.*t);
p_x = -apxPiInvL.*p_X.*sin(apxPiInvL.*x);
p_y = apyPiInvL.*p_Y.*cos(apyPiInvL.*y);
p_z = -apzPiInvL.*p_Z.*sin(apzPiInvL.*z);
p_xx = -apxPiInvL.^2.*p_X.*cos(apxPiInvL.*x);
p_yy = -apyPiInvL.^2.*p_Y.*sin(apyPiInvL.*y);
p_zz = -apzPiInvL.^2.*p_Z.*cos(apzPiInvL.*z);

% Store results in output arrays
q = [r, u, v, w, p];
qt = [r_t, u_t, v_t, w_t, p_t];
qx = [r_x, u_x, v_x, w_x, p_x];
qy = [r_y, u_y, v_y, w_y, p_y];
qz = [r_z, u_z, v_z, w_z, p_z];
qxx = [r_xx, u_xx, u_xy, u_xz, u_yy, u_zz, p_xx];
qyy = [r_yy, v_xx, v_xy, v_yy, v_yz, v_zz, p_yy];
qzz = [r_zz, w_xx, w_xz, w_yy, w_yz, w_zz, p_zz];

%% Source Terms & auxiliary variables
RHO_3 = r;
U_3 = u;
V_3 = v;
W_3 = w;
P_3 = p;

% Source Terms
Q_rho = RHO_3.*auxPiInvL.*u_X.*cos(auxPiInvL.*x) + RHO_3.*avyPiInvL.*v_Y.*cos(avyPiInvL.*y) - RHO_3.*awzPiInvL.*w_Z.*sin(awzPiInvL.*z) + U_3.*arxPiInvL.*r_X.*cos(arxPiInvL.*x) - V_3.*aryPiInvL.*r_Y.*sin(aryPiInvL.*y) + W_3.*arzPiInvL.*r_Z.*cos(arzPiInvL.*z) + artPiInvL.*r_T.*cos(artPiInvL.*t);
Q_rhou = 2*RHO_3.*U_3.*auxPiInvL.*u_X.*cos(auxPiInvL.*x) + RHO_3.*U_3.*avyPiInvL.*v_Y.*cos(avyPiInvL.*y) - RHO_3.*U_3.*awzPiInvL.*w_Z.*sin(awzPiInvL.*z) - RHO_3.*V_3.*auyPiInvL.*u_Y.*sin(auyPiInvL.*y) - RHO_3.*W_3.*auzPiInvL.*u_Z.*sin(auzPiInvL.*z) - RHO_3.*autPiInvL.*u_T.*sin(autPiInvL.*t) + U_3.^2.*arxPiInvL.*r_X.*cos(arxPiInvL.*x) - U_3.*V_3.*aryPiInvL.*r_Y.*sin(aryPiInvL.*y) + U_3.*W_3.*arzPiInvL.*r_Z.*cos(arzPiInvL.*z) + U_3.*artPiInvL.*r_T.*cos(artPiInvL.*t) - apxPiInvL.*p_X.*sin(apxPiInvL.*x) + 1.33333333333333*auxPiInvL.^2.*mu.*u_X.*sin(auxPiInvL.*x) + auyPiInvL.^2.*mu.*u_Y.*cos(auyPiInvL.*y) + auzPiInvL.^2.*mu.*u_Z.*cos(auzPiInvL.*z);
Q_rhov = -RHO_3.*U_3.*avxPiInvL.*v_X.*sin(avxPiInvL.*x) + RHO_3.*V_3.*auxPiInvL.*u_X.*cos(auxPiInvL.*x) + 2*RHO_3.*V_3.*avyPiInvL.*v_Y.*cos(avyPiInvL.*y) - RHO_3.*V_3.*awzPiInvL.*w_Z.*sin(awzPiInvL.*z) + RHO_3.*W_3.*avzPiInvL.*v_Z.*cos(avzPiInvL.*z) + RHO_3.*avtPiInvL.*v_T.*cos(avtPiInvL.*t) + U_3.*V_3.*arxPiInvL.*r_X.*cos(arxPiInvL.*x) - V_3.^2.*aryPiInvL.*r_Y.*sin(aryPiInvL.*y) + V_3.*W_3.*arzPiInvL.*r_Z.*cos(arzPiInvL.*z) + V_3.*artPiInvL.*r_T.*cos(artPiInvL.*t) + apyPiInvL.*p_Y.*cos(apyPiInvL.*y) + avxPiInvL.^2.*mu.*v_X.*cos(avxPiInvL.*x) + 1.33333333333333*avyPiInvL.^2.*mu.*v_Y.*sin(avyPiInvL.*y) + avzPiInvL.^2.*mu.*v_Z.*sin(avzPiInvL.*z);
Q_rhow = RHO_3.*U_3.*awxPiInvL.*w_X.*cos(awxPiInvL.*x) + RHO_3.*V_3.*awyPiInvL.*w_Y.*cos(awyPiInvL.*y) + RHO_3.*W_3.*auxPiInvL.*u_X.*cos(auxPiInvL.*x) + RHO_3.*W_3.*avyPiInvL.*v_Y.*cos(avyPiInvL.*y) - 2*RHO_3.*W_3.*awzPiInvL.*w_Z.*sin(awzPiInvL.*z) - RHO_3.*awtPiInvL.*w_T.*sin(awtPiInvL.*t) + U_3.*W_3.*arxPiInvL.*r_X.*cos(arxPiInvL.*x) - V_3.*W_3.*aryPiInvL.*r_Y.*sin(aryPiInvL.*y) + W_3.^2.*arzPiInvL.*r_Z.*cos(arzPiInvL.*z) + W_3.*artPiInvL.*r_T.*cos(artPiInvL.*t) - apzPiInvL.*p_Z.*sin(apzPiInvL.*z) + awxPiInvL.^2.*mu.*w_X.*sin(awxPiInvL.*x) + awyPiInvL.^2.*mu.*w_Y.*sin(awyPiInvL.*y) + 1.33333333333333*awzPiInvL.^2.*mu.*w_Z.*cos(awzPiInvL.*z);
Q_rhoE = P_3.*auxPiInvL.*u_X.*cos(auxPiInvL.*x) + P_3.*avyPiInvL.*v_Y.*cos(avyPiInvL.*y) - P_3.*awzPiInvL.*w_Z.*sin(awzPiInvL.*z) + RHO_3.*U_3.*(U_3.*auxPiInvL.*u_X.*cos(auxPiInvL.*x) - V_3.*avxPiInvL.*v_X.*sin(avxPiInvL.*x) + W_3.*awxPiInvL.*w_X.*cos(awxPiInvL.*x) + (-P_3.*arxPiInvL.*r_X.*cos(arxPiInvL.*x) - RHO_3.*apxPiInvL.*p_X.*sin(apxPiInvL.*x))./(RHO_3.^2.*(gamma - 1))) + RHO_3.*V_3.*(-U_3.*auyPiInvL.*u_Y.*sin(auyPiInvL.*y) + V_3.*avyPiInvL.*v_Y.*cos(avyPiInvL.*y) + W_3.*awyPiInvL.*w_Y.*cos(awyPiInvL.*y) + (P_3.*aryPiInvL.*r_Y.*sin(aryPiInvL.*y) + RHO_3.*apyPiInvL.*p_Y.*cos(apyPiInvL.*y))./(RHO_3.^2.*(gamma - 1))) + RHO_3.*W_3.*(-U_3.*auzPiInvL.*u_Z.*sin(auzPiInvL.*z) + V_3.*avzPiInvL.*v_Z.*cos(avzPiInvL.*z) - W_3.*awzPiInvL.*w_Z.*sin(awzPiInvL.*z) + (-P_3.*arzPiInvL.*r_Z.*cos(arzPiInvL.*z) - RHO_3.*apzPiInvL.*p_Z.*sin(apzPiInvL.*z))./(RHO_3.^2.*(gamma - 1))) + RHO_3.*auxPiInvL.*u_X.*(P_3./(RHO_3.*(gamma - 1)) + U_3.^2/2 + V_3.^2/2 + W_3.^2/2).*cos(auxPiInvL.*x) + RHO_3.*avyPiInvL.*v_Y.*(P_3./(RHO_3.*(gamma - 1)) + U_3.^2/2 + V_3.^2/2 + W_3.^2/2).*cos(avyPiInvL.*y) - RHO_3.*awzPiInvL.*w_Z.*(P_3./(RHO_3.*(gamma - 1)) + U_3.^2/2 + V_3.^2/2 + W_3.^2/2).*sin(awzPiInvL.*z) + RHO_3.*(-U_3.*autPiInvL.*u_T.*sin(autPiInvL.*t) + V_3.*avtPiInvL.*v_T.*cos(avtPiInvL.*t) - W_3.*awtPiInvL.*w_T.*sin(awtPiInvL.*t) + (-P_3.*artPiInvL.*r_T.*cos(artPiInvL.*t) - RHO_3.*aptPiInvL.*p_T.*sin(aptPiInvL.*t))./(RHO_3.^2.*(gamma - 1))) - U_3.*apxPiInvL.*p_X.*sin(apxPiInvL.*x) + U_3.*arxPiInvL.*r_X.*(P_3./(RHO_3.*(gamma - 1)) + U_3.^2/2 + V_3.^2/2 + W_3.^2/2).*cos(arxPiInvL.*x) + 1.33333333333333*U_3.*auxPiInvL.^2.*mu.*u_X.*sin(auxPiInvL.*x) + U_3.*auyPiInvL.^2.*mu.*u_Y.*cos(auyPiInvL.*y) + U_3.*auzPiInvL.^2.*mu.*u_Z.*cos(auzPiInvL.*z) + V_3.*apyPiInvL.*p_Y.*cos(apyPiInvL.*y) - V_3.*aryPiInvL.*r_Y.*(P_3./(RHO_3.*(gamma - 1)) + U_3.^2/2 + V_3.^2/2 + W_3.^2/2).*sin(aryPiInvL.*y) + V_3.*avxPiInvL.^2.*mu.*v_X.*cos(avxPiInvL.*x) + 1.33333333333333*V_3.*avyPiInvL.^2.*mu.*v_Y.*sin(avyPiInvL.*y) + V_3.*avzPiInvL.^2.*mu.*v_Z.*sin(avzPiInvL.*z) - W_3.*apzPiInvL.*p_Z.*sin(apzPiInvL.*z) + W_3.*arzPiInvL.*r_Z.*(P_3./(RHO_3.*(gamma - 1)) + U_3.^2/2 + V_3.^2/2 + W_3.^2/2).*cos(arzPiInvL.*z) + W_3.*awxPiInvL.^2.*mu.*w_X.*sin(awxPiInvL.*x) + W_3.*awyPiInvL.^2.*mu.*w_Y.*sin(awyPiInvL.*y) + 1.33333333333333*W_3.*awzPiInvL.^2.*mu.*w_Z.*cos(awzPiInvL.*z) + artPiInvL.*r_T.*(P_3./(RHO_3.*(gamma - 1)) + U_3.^2/2 + V_3.^2/2 + W_3.^2/2).*cos(artPiInvL.*t) - auxPiInvL.*u_X.*(2*auxPiInvL.*mu.*u_X.*cos(auxPiInvL.*x) - 0.666666666666667*mu.*(auxPiInvL.*u_X.*cos(auxPiInvL.*x) + avyPiInvL.*v_Y.*cos(avyPiInvL.*y) - awzPiInvL.*w_Z.*sin(awzPiInvL.*z))).*cos(auxPiInvL.*x) + auyPiInvL.*mu.*u_Y.*(-auyPiInvL.*u_Y.*sin(auyPiInvL.*y) - avxPiInvL.*v_X.*sin(avxPiInvL.*x)).*sin(auyPiInvL.*y) + auzPiInvL.*mu.*u_Z.*(-auzPiInvL.*u_Z.*sin(auzPiInvL.*z) + awxPiInvL.*w_X.*cos(awxPiInvL.*x)).*sin(auzPiInvL.*z) + avxPiInvL.*mu.*v_X.*(-auyPiInvL.*u_Y.*sin(auyPiInvL.*y) - avxPiInvL.*v_X.*sin(avxPiInvL.*x)).*sin(avxPiInvL.*x) - avyPiInvL.*v_Y.*(2*avyPiInvL.*mu.*v_Y.*cos(avyPiInvL.*y) - 0.666666666666667*mu.*(auxPiInvL.*u_X.*cos(auxPiInvL.*x) + avyPiInvL.*v_Y.*cos(avyPiInvL.*y) - awzPiInvL.*w_Z.*sin(awzPiInvL.*z))).*cos(avyPiInvL.*y) - avzPiInvL.*mu.*v_Z.*(avzPiInvL.*v_Z.*cos(avzPiInvL.*z) + awyPiInvL.*w_Y.*cos(awyPiInvL.*y)).*cos(avzPiInvL.*z) - awxPiInvL.*mu.*w_X.*(-auzPiInvL.*u_Z.*sin(auzPiInvL.*z) + awxPiInvL.*w_X.*cos(awxPiInvL.*x)).*cos(awxPiInvL.*x) - awyPiInvL.*mu.*w_Y.*(avzPiInvL.*v_Z.*cos(avzPiInvL.*z) + awyPiInvL.*w_Y.*cos(awyPiInvL.*y)).*cos(awyPiInvL.*y) + awzPiInvL.*w_Z.*(-2*awzPiInvL.*mu.*w_Z.*sin(awzPiInvL.*z) - 0.666666666666667*mu.*(auxPiInvL.*u_X.*cos(auxPiInvL.*x) + avyPiInvL.*v_Y.*cos(avyPiInvL.*y) - awzPiInvL.*w_Z.*sin(awzPiInvL.*z))).*sin(awzPiInvL.*z) - gamma.*mu.*(-RHO_3.^2.*apxPiInvL.^2.*p_X.*cos(apxPiInvL.*x) - RHO_3.*(-P_3.*arxPiInvL.^2.*r_X.*sin(arxPiInvL.*x) - 2*apxPiInvL.*arxPiInvL.*p_X.*r_X.*sin(apxPiInvL.*x).*cos(arxPiInvL.*x)) + arxPiInvL.^2.*r_X.^2.*(2*p_0 + 2*p_T.*cos(aptPiInvL.*t) + 2*p_X.*cos(apxPiInvL.*x) + 2*p_Y.*sin(apyPiInvL.*y) + 2*p_Z.*cos(apzPiInvL.*z)).*cos(arxPiInvL.*x).^2)./(Pr.*RHO_3.^3.*(gamma - 1)) - gamma.*mu.*(-RHO_3.^2.*apyPiInvL.^2.*p_Y.*sin(apyPiInvL.*y) - RHO_3.*(-P_3.*aryPiInvL.^2.*r_Y.*cos(aryPiInvL.*y) - 2*apyPiInvL.*aryPiInvL.*p_Y.*r_Y.*sin(aryPiInvL.*y).*cos(apyPiInvL.*y)) + aryPiInvL.^2.*r_Y.^2.*(2*p_0 + 2*p_T.*cos(aptPiInvL.*t) + 2*p_X.*cos(apxPiInvL.*x) + 2*p_Y.*sin(apyPiInvL.*y) + 2*p_Z.*cos(apzPiInvL.*z)).*sin(aryPiInvL.*y).^2)./(Pr.*RHO_3.^3.*(gamma - 1)) - gamma.*mu.*(-RHO_3.^2.*apzPiInvL.^2.*p_Z.*cos(apzPiInvL.*z) - RHO_3.*(-P_3.*arzPiInvL.^2.*r_Z.*sin(arzPiInvL.*z) - 2*apzPiInvL.*arzPiInvL.*p_Z.*r_Z.*sin(apzPiInvL.*z).*cos(arzPiInvL.*z)) + arzPiInvL.^2.*r_Z.^2.*(2*p_0 + 2*p_T.*cos(aptPiInvL.*t) + 2*p_X.*cos(apxPiInvL.*x) + 2*p_Y.*sin(apyPiInvL.*y) + 2*p_Z.*cos(apzPiInvL.*z)).*cos(arzPiInvL.*z).^2)./(Pr.*RHO_3.^3.*(gamma - 1));

source = [Q_rho, Q_rhou, Q_rhov, Q_rhow, Q_rhoE];
end % funtion