clc
clear

% draw the manifold 


% x = -10:1:10;
% y = -10:1:10;
% 
% pphipx = [1,1];
% 
% z = [];
% for i = 1:size(x,2)
%     xx = x(i);
%     yy = y(i);
%     zz = pphipx * [xx,yy]';
%     z = [z zz];
% end

[X,Y] = meshgrid(-10:1:10,-10:1:10);
Z = X+Y;

surf(X,Y,Z)