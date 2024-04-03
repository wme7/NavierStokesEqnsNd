function flag = evalNavierStokesEqns3d(q, qt, qx, qy, qz, qxx, qyy, qzz, source, TOL)

% Flow parameters
global gamma mu Pr %#ok<GVMIS>

% Solutions and derivatives
r=q(:,1); r_t=qt(:,1); r_x=qx(:,1); r_y=qy(:,1); r_z=qz(:,1); r_xx=qxx(:,1);                               r_yy=qyy(:,1); r_zz=qzz(:,1);
u=q(:,2); u_t=qt(:,2); u_x=qx(:,2); u_y=qy(:,2); u_z=qz(:,2); u_xx=qxx(:,2); u_xy=qxx(:,3); u_xz=qxx(:,4); u_yy=qxx(:,5); u_zz=qxx(:,6);
v=q(:,3); v_t=qt(:,3); v_x=qx(:,3); v_y=qy(:,3); v_z=qz(:,3); v_xx=qyy(:,2); v_xy=qyy(:,3); v_yz=qyy(:,4); v_yy=qyy(:,5); v_zz=qyy(:,6);
w=q(:,4); w_t=qt(:,4); w_x=qx(:,4); w_y=qy(:,4); w_z=qz(:,4); w_xx=qzz(:,2); w_xz=qzz(:,3); w_yz=qzz(:,4); w_yy=qzz(:,5); w_zz=qzz(:,6);
p=q(:,5); p_t=qt(:,5); p_x=qx(:,5); p_y=qy(:,5); p_z=qz(:,5); p_xx=qxx(:,7);                               p_yy=qyy(:,7); p_zz=qzz(:,7);

e = 1 ./ (gamma - 1) .* (p ./ r);
e_t = 1 ./ (gamma - 1) .* (p_t .* r - p .* r_t) ./ r.^2;
e_x = 1 ./ (gamma - 1) .* (p_x .* r - p .* r_x) ./ r.^2;
e_y = 1 ./ (gamma - 1) .* (p_y .* r - p .* r_y) ./ r.^2;
e_z = 1 ./ (gamma - 1) .* (p_z .* r - p .* r_z) ./ r.^2;
e_xx = 1 ./ (gamma - 1) .* (p_xx .* r.^2 - r .* (2 .* p_x .* r_x + p .* r_xx) + 2 .* p .* r_x.^2) ./ r.^3;
e_yy = 1 ./ (gamma - 1) .* (p_yy .* r.^2 - r .* (2 .* p_y .* r_y + p .* r_yy) + 2 .* p .* r_y.^2) ./ r.^3;
e_zz = 1 ./ (gamma - 1) .* (p_zz .* r.^2 - r .* (2 .* p_z .* r_z + p .* r_zz) + 2 .* p .* r_z.^2) ./ r.^3;

E = e + (u.^2 + v.^2 + w.^2) / 2;
E_t = e_t + u .* u_t + v .* v_t + w .* w_t;
E_x = e_x + u .* u_x + v .* v_x + w .* w_x;
E_y = e_y + u .* u_y + v .* v_y + w .* w_y;
E_z = e_z + u .* u_z + v .* v_z + w .* w_z;

% mu = constant
mu_x = 0;
mu_y = 0;
mu_z = 0;

% Assuming Stokes' Hypothesis
lambda = -2/3 * mu;
lambda_x = 0;
lambda_y = 0;
lambda_z = 0;

% Define heat flux
% q_x = -(gamma / Pr) .* mu .* e_x;
% q_y = -(gamma / Pr) .* mu .* e_y;
% q_z = -(gamma / Pr) .* mu .* e_z;
q_xx = -(gamma / Pr) .* (mu_x .* e_x + mu .* e_xx);
q_yy = -(gamma / Pr) .* (mu_y .* e_y + mu .* e_yy);
q_zz = -(gamma / Pr) .* (mu_z .* e_z + mu .* e_zz);

% Compute some quantities needed to evaluate the 2D unsteady Euler system
ru_t = r_t .* u + r .* u_t;
rv_t = r_t .* v + r .* v_t;
rw_t = r_t .* w + r .* w_t;
ru_x = r_x .* u + r .* u_x;
rv_y = r_y .* v + r .* v_y;
rw_z = r_z .* w + r .* w_z;
rE_t = r_t .* E + r .* E_t;
ruu_x = (r_x .* u .* u) + (r .* u_x .* u) + (r .* u .* u_x);
ruv_x = (r_x .* u .* v) + (r .* u_x .* v) + (r .* u .* v_x);
ruw_x = (r_x .* w .* u) + (r .* w_x .* u) + (r .* w .* u_x);
ruv_y = (r_y .* u .* v) + (r .* u_y .* v) + (r .* u .* v_y);
rvv_y = (r_y .* v .* v) + (r .* v_y .* v) + (r .* v .* v_y);
rwv_y = (r_y .* w .* v) + (r .* w_y .* v) + (r .* w .* v_y);
ruw_z = (r_z .* u .* w) + (r .* u_z .* w) + (r .* u .* w_z);
rvw_z = (r_z .* v .* w) + (r .* v_z .* w) + (r .* v .* w_z);
rww_z = (r_z .* w .* w) + (r .* w_z .* w) + (r .* w .* w_z);
ruE_x = (r_x .* u .* E) + (r .* u_x .* E) + (r .* u .* E_x);
rvE_y = (r_y .* v .* E) + (r .* v_y .* E) + (r .* v .* E_y);
rwE_z = (r_z .* w .* E) + (r .* w_z .* E) + (r .* w .* E_z);

% Define stress tensor variables
tauxx = 2 .* mu .* u_x + lambda .* (u_x + v_y + w_z);
tauyy = 2 .* mu .* v_y + lambda .* (u_x + v_y + w_z);
tauzz = 2 .* mu .* w_z + lambda .* (u_x + v_y + w_z);
tauxy = mu .* (u_y + v_x);
tauxz = mu .* (u_z + w_x);
tauyz = mu .* (v_z + w_y);

tauxx_x = 2 .* mu_x .* u_x + 2 .* mu .* u_xx + lambda_x .* (u_x + v_y + w_z) + lambda .* (u_xx + v_xy + w_xz);
tauyy_y = 2 .* mu_y .* v_y + 2 .* mu .* v_yy + lambda_y .* (u_x + v_y + w_z) + lambda .* (u_xy + v_yy + w_yz);
tauzz_z = 2 .* mu_z .* w_z + 2 .* mu .* w_zz + lambda_z .* (u_x + v_y + w_z) + lambda .* (u_xz + v_yz + w_zz);

tauxy_x = mu_x .* (u_y + v_x) + mu .* (u_xy + v_xx);
tauxy_y = mu_y .* (u_y + v_x) + mu .* (u_yy + v_xy);
tauxz_x = mu_x .* (u_z + w_x) + mu .* (u_xz + w_xx);
tauxz_z = mu_z .* (u_z + w_x) + mu .* (u_zz + w_xz);
tauyz_y = mu_y .* (v_z + w_y) + mu .* (v_yz + w_yy);
tauyz_z = mu_z .* (v_z + w_y) + mu .* (v_zz + w_yz);

% Pressure x velocities
pu_x = p_x .* u + p .* u_x;
pv_y = p_y .* v + p .* v_y;
pw_z = p_z .* w + p .* w_z;

% Tau components x velocities
utauxx_x = u_x .* tauxx + u .* tauxx_x;
vtauxy_x = v_x .* tauxy + v .* tauxy_x;
wtauxz_x = w_x .* tauxz + w .* tauxz_x;
utauxy_y = u_y .* tauxy + u .* tauxy_y;
vtauyy_y = v_y .* tauyy + v .* tauyy_y;
wtauyz_y = w_y .* tauyz + w .* tauyz_y;
utauxz_z = u_z .* tauxz + u .* tauxz_z;
vtauyz_z = v_z .* tauyz + v .* tauyz_z;
wtauzz_z = w_z .* tauzz + w .* tauzz_z;

% Form the 2D unsteady Euler equations
equation(1) = sum(  r_t +  ru_x +  rv_y +  rw_z  - source(:,1));
equation(2) = sum( ru_t + ruu_x + ruv_y + ruw_z + p_x - (tauxx_x + tauxy_y + tauxz_z) - source(:,2));
equation(3) = sum( rv_t + ruv_x + rvv_y + rvw_z + p_y - (tauxy_x + tauyy_y + tauyz_z) - source(:,3));
equation(4) = sum( rw_t + ruw_x + rwv_y + rww_z + p_z - (tauxz_x + tauyz_y + tauzz_z) - source(:,4));
equation(5) = sum( rE_t + ruE_x + rvE_y + rwE_z + pu_x + pv_y + pw_z ...
    - utauxx_x - vtauxy_x - wtauxz_x ...
    - utauxy_y - vtauyy_y - wtauyz_y ...
    - utauxz_z - vtauyz_z - wtauzz_z ...
    + q_xx + q_yy + q_zz - source(:,5));

% Display the results. All must be zero.
fprintf('\nSubstitution yields:\n\n');
fprintf('Continuity = %1.12f\n', equation(1));
fprintf('X-momentum = %1.12f\n', equation(2));
fprintf('Y-momentum = %1.12f\n', equation(3));
fprintf('Z-momentum = %1.12f\n', equation(4));
fprintf('Energy     = %1.12f\n', equation(5));

% Output
flag = all(equation < TOL);
if (flag)
    fprintf('\n Is an exact solution!\n\n'); 
else
    fprintf('\n Not an exact solution!\n\n'); 
end

end % function