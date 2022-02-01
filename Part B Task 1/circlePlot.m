function circlePlot(pos, rho)
% pos (3x1 vector) = position of Rx
% rho (float) = distance betwwen one Rx and Tx

temp = 0:pi/100:2*pi;
x = rho * cos(temp) + pos(1);
y = rho * sin(temp) + pos(2);
plot(x, y, pos(1), pos(2),'ok','MarkerFaceColor','y');
end