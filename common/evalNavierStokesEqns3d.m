function flag = evalNavierStokesEqns3d(q, qx, qy, qz, qxx, qyy, qzz, qt, TOL)

% Flow parameters
global gamma mu Pr

% Solutions and derivatives
r = q(:,1); invr = 1./r;
u = q(:,2);
v = q(:,3);
w = q(:,4);
p = q(:,5);

rx = qx(:,1);
ux = qx(:,2);
vx = qx(:,3);
wx = qx(:,4);
px = qx(:,5);

ry = qy(:,1);
uy = qy(:,2);
vy = qy(:,3);
wy = qy(:,4);
py = qy(:,5);

rz = qz(:,1);
uz = qz(:,2);
vz = qz(:,3);
wz = qz(:,4);
pz = qz(:,5);

rxx = qxx(:,1);
uxx = qxx(:,2);
uxy = qxx(:,3);
uxz = qxx(:,4);
uyy = qxx(:,5);
uzz = qxx(:,6);
pxx = qxx(:,7);

ryy = qyy(:,1);
vxx = qyy(:,2);
vxy = qyy(:,3);
vyy = qyy(:,4);
vyz = qyy(:,5);
vzz = qyy(:,6);
pyy = qyy(:,7);

rzz = qzz(:,1);
wxx = qzz(:,2);
wyy = qzz(:,3);
wxz = qzz(:,4);
wyz = qzz(:,5);
wzz = qzz(:,6);
pzz = qzz(:,7);

rt = qt(:,1);
ut = qt(:,2);
vt = qt(:,3);
wt = qt(:,4);
pt = qt(:,5);

% Compute some quantities needed to evaluate the 2D unsteady Euler system
rhout = rt .* u + r .* ut;
rhovt = rt .* v + r .* vt;
rhowt = rt .* w + r .* wt;

rhoux = rx .* u + r .* ux;
rhovy = ry .* v + r .* vy;
rhowz = rz .* w + r .* wz;

rhouux = rx .* u .* u + r .* ux .* u + r .* u .* ux;
rhouvx = rx .* v .* u + r .* vx .* u + r .* v .* ux;
rhouwx = rx .* w .* u + r .* wx .* u + r .* w .* ux;

rhouvy = ry .* u .* v + r .* uy .* v + r .* u .* vy;
rhovvy = ry .* v .* v + r .* vy .* v + r .* v .* vy;
rhowvy = ry .* w .* v + r .* wy .* v + r .* w .* vy;

rhouwz = rz .* u .* w + r .* uz .* w + r .* u .* wz;
rhovwz = rz .* v .* w + r .* vz .* w + r .* v .* wz;
rhowwz = rz .* w .* w + r .* wz .* w + r .* w .* wz;

rH  = gamma * p  ./ (gamma - 1) + 0.5 *   r .* (u .* u + v .* v + w .* w);
rHx = gamma * px ./ (gamma - 1) + 0.5 * (rx .* (u .* u + v .* v + w .* w) + 2 * r .* (u .* ux + v .* vx + w .* wx));
rHy = gamma * py ./ (gamma - 1) + 0.5 * (ry .* (u .* u + v .* v + w .* w) + 2 * r .* (u .* uy + v .* vy + w .* wy));
rHz = gamma * pz ./ (gamma - 1) + 0.5 * (rz .* (u .* u + v .* v + w .* w) + 2 * r .* (u .* uz + v .* vz + w .* wz));

rhouHx = ux .* rH + u .* rHx;
rhovHy = vy .* rH + v .* rHy;
rhowHz = wz .* rH + w .* rHz;

rEt = pt ./ (gamma - 1) + 0.5 * (rt .* (u .* u + v .* v + w .* w) + 2 * r .* (u .* ut + v .* vt + w .* wt));

% Stokes' hypothesis
lmbd = -2/3 * mu;

% Viscous Stress tensor
tau11 = 2 * mu * ux + lmbd * (ux + vy + wz);
tau22 = 2 * mu * vy + lmbd * (ux + vy + wz);
tau33 = 2 * mu * wz + lmbd * (ux + vy + wz);
tau12 = mu * (uy + vx);       tau21 = tau12;
tau13 = mu * (wx + uz);       tau31 = tau13;
tau23 = mu * (vz + wy);       tau32 = tau23;
tau11x = 2 * mu * uxx + lmbd * (uxx + vxy + wxz);
tau22y = 2 * mu * vyy + lmbd * (uxy + vyy + wyz);
tau33z = 2 * mu * wzz + lmbd * (uxz + vyz + wzz);
tau12x = mu * (uxy + vxx); %tau21x = tau12x;
tau12y = mu * (uyy + vxy); tau21y = tau12y;
tau13x = mu * (wxx + uxz); %tau31x = tau13x;
tau13z = mu * (wxz + uzz); tau31z = tau13z;
tau23y = mu * (vyz + wyy); %tau32y = tau23y;
tau23z = mu * (vzz + wyz); tau32z = tau23z;

% Heat Fluxes
q1x = gamma * mu /(Pr * (gamma-1)) * (-2.*invr.*invr.* rx.* px + invr.* pxx + p.*invr.*invr.*invr.* rx.* rx - p.*invr.*invr.* rxx);
q2y = gamma * mu /(Pr * (gamma-1)) * (-2.*invr.*invr.* ry.* py + invr.* pyy + p.*invr.*invr.*invr.* ry.* ry - p.*invr.*invr.* ryy);
q3z = gamma * mu /(Pr * (gamma-1)) * (-2.*invr.*invr.* rz.* pz + invr.* pzz + p.*invr.*invr.*invr.* rz.* rz - p.*invr.*invr.* rzz);

tauU1x = u.* tau11x + v.* tau12x + w.* tau13x + tau11.* ux + tau12.* vx + tau13.* wx;
tauU2y = u.* tau21y + v.* tau22y + w.* tau23y + tau21.* uy + tau22.* vy + tau23.* wy;
tauU3z = u.* tau31z + v.* tau32z + w.* tau33z + tau31.* uz + tau32.* vz + tau33.* wz;

% Form the 2D unsteady Euler equations
equation(1) = sum(    rt +  rhoux +  rhovy + rhowz );
equation(2) = sum( rhout + rhouux + rhouvy + rhouwz + px - tau11x - tau12x - tau13x);
equation(3) = sum( rhovt + rhouvx + rhovvy + rhovwz + py - tau21y - tau22y - tau23y);
equation(4) = sum( rhowt + rhouwx + rhowvy + rhowwz + pz - tau31z - tau32z - tau33z);
equation(5) = sum(   rEt + rhouHx + rhovHy + rhowHz - tauU1x - tauU2y - tauU3z + q1x + q2y + q3z);

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