function [hr, hu, hv, hw, hp] = plotNavierStokesEqns3d(x,y,z,q,t)
    DIM = 3;
    colormap jet;
    L = max(x(:));
    subplot(231);
        hr = slice(x,y,z,reshape(q(:,1),size(z)),[0,L],[0,L],[0,L]);
        view(DIM);
        title(sprintf('$\\rho(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        zlabel('$z$', Interpreter='latex');
        shading interp;
        clim([0,2]);
        colorbar;
    subplot(232);
        hu = slice(x,y,z,reshape(q(:,2),size(z)),[0,L],[0,L],[0,L]);
        view(DIM);
        title(sprintf('$u(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        zlabel('$z$', Interpreter='latex');
        shading interp;
        clim([-1,1]);
        colorbar;
    subplot(233);
        hv = slice(x,y,z,reshape(q(:,3),size(z)),[0,L],[0,L],[0,L]);
        view(DIM);
        title(sprintf('$v(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        zlabel('$z$', Interpreter='latex');
        shading interp;
        clim([-1,1]);
        colorbar;
    subplot(234);
        hw = slice(x,y,z,reshape(q(:,4),size(z)),[0,L],[0,L],[0,L]);
        view(DIM);
        title(sprintf('$w(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        zlabel('$z$', Interpreter='latex');
        shading interp;
        clim([-1,1]);
        colorbar;
    subplot(235);
        hp = slice(x,y,z,reshape(q(:,5),size(z)),[0,L],[0,L],[0,L]);
        view(DIM);
        title(sprintf('$\\wp(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        zlabel('$z$', Interpreter='latex');
        shading interp;
        clim([0,2]);
        colorbar;
end
