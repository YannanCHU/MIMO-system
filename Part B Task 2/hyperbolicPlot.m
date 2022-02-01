function hyperbolicPlot(r1, ri, rhoi1)
% r1 (3x1 vector) = position of Rx 1
% ri (3x1 vector) = position of Rx i
% rhoi1 (float) = difference between the length from r1 to Tx and the
% length from ri to Tx

f = @(x,y) abs(sqrt((x-r1(1)).^2 + (y-r1(2)).^2) - sqrt((x-ri(1)).^2 + (y-ri(2)).^2)).* ...
    (sqrt((x-r1(1)).^2 + (y-r1(2)).^2) - sqrt((x-ri(1)).^2 + (y-ri(2)).^2) < 0) - ...
    rhoi1;
fimplicit(f);
plot(ri(1), ri(2),'ok','MarkerFaceColor','y');
end