% Plot the final experimental designs

grid_interv = 0.01;
[xs1, xs2] = meshgrid(-8:grid_interv:8);
G_true = G_function(xs1, xs2); % True function on the grid points
[mu_grid, s2_grid] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x_doe, G_doe, [xs1(:) xs2(:)]);
mu_grid = reshape(mu_grid, size(G_true)); s2_grid = reshape(s2_grid, size(G_true));

radi_dp = 3;
theta_circ = linspace(0,2*pi);
x_circleplot = radi_dp*cos(theta_circ);
y_circleplot = radi_dp*sin(theta_circ);
radi_dp2 = 3.5;
x_circleplot2 = radi_dp2*cos(theta_circ);
y_circleplot2= radi_dp2*sin(theta_circ);

Finfig = figure();
Finfig.Color = 'w';
Finfig.Position = [680 430 760 475];
contour(xs1, xs2, G_true,'levellist',0,'LineStyle','-', 'LineWidth', 2.0, 'Color',[0.39, 0.47, 0.64]); hold on; % True
contour(xs1, xs2, mu_grid,'levellist',0,'LineStyle','--','LineWidth', 1.5,'Color','k'); % Kriging
plot(x_circleplot,y_circleplot,'LineStyle','-.', 'LineWidth', 1,'Color', [0.40 0.40 0.40],'HandleVisibility','off');
plot(x_circleplot2,y_circleplot2,'LineStyle','-.', 'LineWidth', 1,'Color', [0.40 0.40 0.40],'HandleVisibility','off');
scatter(x_doe(1:N_iniDOE,1),x_doe(1:N_iniDOE,2),30,'k','LineWidth',1.1)
scatter(x_doe(N_iniDOE+1:end,1),x_doe(N_iniDOE+1:end,2),50,'r','x','LineWidth',1.5)
xlabel('x_1','fontsize',11)
ylabel('x_2','fontsize',11)
legend('True function','Kriging prediction','Initial DoE','Added DoE','Location','NorthEast')
axis equal; axis([-7 7 -7 7]); hold off;


