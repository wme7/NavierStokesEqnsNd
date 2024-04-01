function flag = evalNavierStokesEqns2d(q, qt, qx, qy, qxx, qyy, source, TOL)

% Flow parameters
global gamma mu Pr %#ok<GVMIS>

% Solutions and derivatives
r=q(:,1); r_t=qt(:,1); r_x=qx(:,1); r_y=qy(:,1); r_xx=qxx(:,1);                r_yy=qyy(:,1);
u=q(:,2); u_t=qt(:,2); u_x=qx(:,2); u_y=qy(:,2); u_xx=qxx(:,2); u_xy=qxx(:,3); u_yy=qxx(:,4);
v=q(:,3); v_t=qt(:,3); v_x=qx(:,3); v_y=qy(:,3); v_xx=qyy(:,2); v_xy=qyy(:,3); v_yy=qyy(:,4);
p=q(:,4); p_t=qt(:,4); p_x=qx(:,4); p_y=qy(:,4); p_xx=qxx(:,5);                p_yy=qyy(:,5);

% Auxiliary variables
e = 1 / (gamma - 1) * (p ./ r);
e_t = 1 / (gamma - 1) * (p_t .* r - p .* r_t) ./ r.^2;
e_x = 1 / (gamma - 1) * (p_x .* r - p .* r_x) ./ r.^2;
e_y = 1 / (gamma - 1) * (p_y .* r - p .* r_y) ./ r.^2;
e_xx = 1 / (gamma - 1) * (p_xx .* r.^2 - r .* (2 .* p_x .* r_x + p .* r_xx) + 2 .* p .* r_x.^2) ./ r.^3;
e_yy = 1 / (gamma - 1) * (p_yy .* r.^2 - r .* (2 .* p_y .* r_y + p .* r_yy) + 2 .* p .* r_y.^2) ./ r.^3;

E = e + (u.^2 + v.^2) / 2;
E_t = e_t + u .* u_t + v .* v_t;
E_x = e_x + u .* u_x + v .* v_x;
E_y = e_y + u .* u_y + v .* v_y;

% mu = constant
mu_x = 0;
mu_y = 0;

% Assuming Stokes' Hypothesis
lambda = -2/3 * mu;
lambda_x = 0;
lambda_y = 0;

% Define heat flux
% q_x = -(gamma / Pr) .* mu .* e_x;
% q_y = -(gamma / Pr) .* mu .* e_y;
q_xx = -(gamma / Pr) .* (mu_x .* e_x + mu .* e_xx);
q_yy = -(gamma / Pr) .* (mu_y .* e_y + mu .* e_yy);

% Compute some quantities needed to evaluate the 2D unsteady Euler system
ru_t = r_t .* u + r .* u_t;
rv_t = r_t .* v + r .* v_t;
ru_x = r_x .* u + r .* u_x;
rv_y = r_y .* v + r .* v_y;
rE_t = r_t .* E + r .* E_t;
ruu_x = (r_x .* u .* u) + (r .* u_x .* u) + (r .* u .* u_x);
ruv_x = (r_x .* u .* v) + (r .* u_x .* v) + (r .* u .* v_x);
ruv_y = (r_y .* u .* v) + (r .* u_y .* v) + (r .* u .* v_y);
rvv_y = (r_y .* v .* v) + (r .* v_y .* v) + (r .* v .* v_y);
ruE_x = (r_x .* u .* E) + (r .* u_x .* E) + (r .* u .* E_x);
rvE_y = (r_y .* v .* E) + (r .* v_y .* E) + (r .* v .* E_y);

% Define stress tensor variables
tauxx = 2 * mu .* u_x + lambda .* (u_x + v_y);
tauyy = 2 * mu .* v_y + lambda .* (u_x + v_y);
tauxy = mu .* (u_y + v_x);

tauxx_x = 2 .* mu_x .* u_x + 2 .* mu .* u_xx + lambda_x .* (u_x + v_y) + lambda .* (u_xx + v_xy);
tauyy_y = 2 .* mu_y .* v_y + 2 .* mu .* v_yy + lambda_y .* (u_x + v_y) + lambda .* (u_xy + v_yy);

tauxy_x = mu_x .* (u_y + v_x) + mu .* (u_xy + v_xx);
tauxy_y = mu_y .* (u_y + v_x) + mu .* (u_yy + v_xy);

% Pressure x velocities
pu_x = p_x .* u + p .* u_x;
pv_y = p_y .* v + p .* v_y;

% Tau components x velocities
utauxx_x = u_x .* tauxx + u .* tauxx_x;
vtauxy_x = v_x .* tauxy + v .* tauxy_x;
utauxy_y = u_y .* tauxy + u .* tauxy_y;
vtauyy_y = v_y .* tauyy + v .* tauyy_y;

% Form the 2D unsteady Euler equations
equation(1) = sum(  r_t +  ru_x + rv_y - source(:,1));
equation(2) = sum( ru_t + ruu_x + ruv_y + p_x - (tauxx_x + tauxy_y) - source(:,2));
equation(3) = sum( rv_t + ruv_x + rvv_y + p_y - (tauxy_x + tauyy_y) - source(:,3));
equation(4) = sum( rE_t + ruE_x + rvE_y + pu_x + pv_y ...
    - utauxx_x - vtauxy_x ...
    - utauxy_y - vtauyy_y ...
    + q_xx + q_yy - source(:,4));

% Display the results. All must be zero.
fprintf('\nSubstitution yields:\n\n');
fprintf('Continuity = %1.12f\n', equation(1));
fprintf('X-momentum = %1.12f\n', equation(2));
fprintf('Y-momentum = %1.12f\n', equation(3));
fprintf('Energy     = %1.12f\n', equation(4));

% Output
flag = any(equation < TOL);
if (flag)
    fprintf('\n Is an exact solution!\n\n'); 
else
    fprintf('\n Not an exact solution!\n\n'); 
end

end % function