% function f to interpolate
f = @(x,y) (sin(10*x)+cos(y*y*3))*cos(40*y)*(1-exp(x^4));

% grid1 with fineness level 8
grid1 = sasg(2,8);
grid1.hierarchize(f);

% grid2 with fineness level 8
grid2 = sasg(2,8);
grid2.hierarchize(f);
% adaptively refine grid2
grid2.unrerefine(f, 1e-2);

% function handles for plots
g1 = @(x,y) grid1.eval(x,y);
g2 = @(x,y) grid2.eval(x,y);

% plot function f
figure; subplot(2,3,1)
fsurf(f,[0 1 0 1])
zlim([-2 2]); caxis([-2 2])
% plot function g1
subplot(2,3,2)
fsurf(g1,[0 1 0 1])
zlim([-2 2]); caxis([-2 2])
% plot function g2
subplot(2,3,3)
fsurf(g2,[0 1 0 1])
zlim([-2 2]); caxis([-2 2])
% plot grid1 nodes
subplot(2,3,5)
nodes1 = grid1.listnodes();
scatter(nodes1(1,:),nodes1(2,:),'.');
% plot grid2 nodes
subplot(2,3,6)
nodes2 = grid2.listnodes();
scatter(nodes2(1,:),nodes2(2,:),'.');