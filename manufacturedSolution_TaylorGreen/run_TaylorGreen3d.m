% Test Exact solutions
addpath('../common/');
global gamma mu Pr

% Aceptable numerical error to stablish solution validity:
TOL = 1E-12;

% Domain length
L = 1.0;

% Fluid parameters (NOTE: uniform viscosity assumed!)
gamma = 1.4;
mu = 0.0025;
Pr = 0.72;

% 3D - mesh
[x,y,z] = meshgrid(linspace(0,L,45));

% Initialize figure
figure(2)

% Compute the exact solution and derivatives
t0=0; dt=0.05; tEnd=0.5;
for t = t0:dt:tEnd

    % Compute exact solution
    [q,qx,qy,qz,qxx,qyy,qzz,qt] = manufacturedSolution_TaylorGreen3d(x(:),y(:),z(:),L,t);

    % Stop if solution does not meet standard
    if not(evalNavierStokesEqns3d(q,qx,qy,qz,qxx,qyy,qzz,qt,TOL)), return; end

    % Update figure
    % plotNavierStokesEqns3d(x,y,z,q,t);
    plotVorticityNavierStokesEqns3d(x,y,z,q,t)
    drawnow;
end
