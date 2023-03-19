%% eval the maximum sinx*cosx 
clc 
clear
vals = [];
for theta = pi/6:0.001:(5/6*pi)
    vals = [vals sin(theta)*cos(theta)];
end

plot(vals,'o',LineWidth=1.5);
fprintf('the maximum val is %d', max(vals));