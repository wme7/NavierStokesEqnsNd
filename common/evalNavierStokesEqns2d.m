function flag = evalNavierStokesEqns2d(q, qx, qy, qxx, qyy, qt, TOL)

% Flow parameters
global gamma mu Pr

% Solutions and derivatives
r = q(:,1); invr = 1./r;
u = q(:,2);
v = q(:,3);
p = q(:,4);

rx = qx(:,1);
ux = qx(:,2);
vx = qx(:,3);
px = qx(:,4);

ry = qy(:,1);
uy = qy(:,2);
vy = qy(:,3);
py = qy(:,4);

rxx = qxx(:,1);
uxx = qxx(:,2);
uxy = qxx(:,3);
uyy = qxx(:,4);
pxx = qxx(:,5);

ryy = qyy(:,1);
vxx = qyy(:,2);
vxy = qyy(:,3);
vyy = qyy(:,4);
pyy = qyy(:,5);

rt = qt(:,1);
ut = qt(:,2);
vt = qt(:,3);
pt = qt(:,4);

% Compute some quantities needed to evaluate the 2D unsteady Euler system
rhout = rt .* u + r .* ut;
rhovt = rt .* v + r .* vt;

rhoux = rx .* u + r .* ux;
rhovy = ry .* v + r .* vy;

rhouux = rx .* u .* u + r .* ux .* u + r .* u .* ux;
rhouvx = rx .* u .* v + r .* ux .* v + r .* u .* vx;

rhouvy = ry .* u .* v + r .* uy .* v + r .* u .* vy;
rhovvy = ry .* v .* v + r .* vy .* v + r .* v .* vy;

rH  = gamma * p  ./ (gamma - 1) + 0.5 *  r .* (u .* u + v .* v);
rHx = gamma * px ./ (gamma - 1) + 0.5 * (rx .* (u .* u + v .* v) + 2 * r .* (u .* ux + v .* vx));
rHy = gamma * py ./ (gamma - 1) + 0.5 * (ry .* (u .* u + v .* v) + 2 * r .* (u .* uy + v .* vy));

rhouHx = ux .* rH + u .* rHx;
rhovHy = vy .* rH + v .* rHy;

rEt = pt ./ (gamma - 1) + 0.5 * (rt .* (u .* u + v .* v) + 2 * r .* (u .* ut + v .* vt));

% Stokes' hypothesis
lmbd = -2/3 * mu;

% Viscous Stress tensor
tau11 = 2 * mu * ux + lmbd * (ux + vy);
tau22 = 2 * mu * vy + lmbd * (ux + vy);
tau12 = mu * (uy + vx);  tau21 = tau12;
tau11x = 2 * mu * uxx + lmbd * (uxx + vxy);
tau22y = 2 * mu * vyy + lmbd * (uxy + vyy);
tau12x = mu * (uxy + vxx); %tau21x = tau12x;
tau12y = mu * (uyy + vxy); tau21y = tau12y;

% Heat Fluxes
q1x = gamma * mu /(Pr * (gamma-1)) * (-2*invr.*invr.* rx.* px + invr.* pxx + p.*invr.*invr.*invr.* rx.* rx - p.*invr.*invr.* rxx);
q2y = gamma * mu /(Pr * (gamma-1)) * (-2*invr.*invr.* ry.* py + invr.* pyy + p.*invr.*invr.*invr.* ry.* ry - p.*invr.*invr.* ryy);

tauU1x = u.* tau11x + v.* tau12x + tau11.* ux + tau12.* vx;
tauU2y = u.* tau21y + v.* tau22y + tau21.* uy + tau22.* vy;

% Form the 2D unsteady Euler equations
equation(1) = sum(    rt +  rhoux + rhovy );
equation(2) = sum( rhout + rhouux + rhouvy + px - tau11x - tau12x);
equation(3) = sum( rhovt + rhouvx + rhovvy + py - tau21y - tau22y);
equation(4) = sum(   rEt + rhouHx + rhovHy - tauU1x - tauU2y + q1x + q2y);

% Display the results. All must be zero.
fprintf('\nSubstitution yields:\n\n');
fprintf('Continuity = %1.12f\n', equation(1));
fprintf('X-momentum = %1.12f\n', equation(2));
fprintf('Y-momentum = %1.12f\n', equation(3));
fprintf('Energy     = %1.12f\n', equation(4));

% Output
flag = all(equation < TOL);
if (flag)
    fprintf('\n Is an exact solution!\n\n'); 
else
    fprintf('\n Not an exact solution!\n\n'); 
end

end % function