% Test Exact solutions
addpath('../common/');

% Aceptable numerical error to stablish solution validity:
TOL = 1E-12;

% Domain length
L = 1.0;

% 2D - mesh
[x,y] = meshgrid(linspace(0,L,25));

% Initialize figure
figure(1)

% Compute the exact solution and derivatives
t0=0; dt=0.05; tEnd=1;
for t = t0:dt:tEnd
    
    % Compute exact solution
    [q,qx,qy,qt] = manufacturedSolution_TaylorGreen2d(x(:),y(:),L,t);
    
    % Stop if solution does not meet standard
    % if not(evalNavierStokesEqns2d(q,qx,qy,qt,TOL)), return; end

    % Update figure
    % plotNavierStokesEqns2d(x,y,q,t);
    plotVorticityNavierStokesEqns2d(x,y,q,t);
    drawnow;
end