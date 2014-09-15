x = 1:0.1:10;
y1 = sin(x);
y2 = cos(x);
y3 = 100*sin(2*x);
y4 = 100*cos(2*x);

[ax,h1,h2] = plotyy(x,y1,x,y3);

% Add additional plot to the respective axes
hold(ax(1));
plot(ax(1),x,y2)
hold(ax(2));
plot(ax(2),x,y4)

% Creates a common legend
legend('y1','y3','y2','y4')

% % Creates two seperate legends
% legend(ax(1),'y1','y3','Location','NorthWest')
% legend(ax(2),'y2','y4','Location','NorthEast')