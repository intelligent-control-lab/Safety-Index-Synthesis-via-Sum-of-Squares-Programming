clc
clear 

%% initialization  
max_iter = 1000000;
x = [pi/2,0]; % initial state, theta = pi/2, dtheta = 0.
dt = 1e-6;
% k = 3.741078e-06;
% k = 0.5324;
% k = 10;
% k = 3.3634;
% k = 178;
% k = 2.63;
k = 2
save_freq = 20;

thetalim = [0, pi];
dthetalim = [-1,1];
ddthetalim = [-1, 1];


%% simulation 
phis = [];
for i = 1:max_iter
%     u_r = 2*rand(1) - 1; % random reference control
    u_r = -1; % the maximum unsafe control 
    theta = x(1);
    dtheta = x(2);
    if mod(i, save_freq) == 0
        phis = [phis phi(x)];
    end
    if phi(x) > -1e-5 % due to discrete time system 
        % apply the maximum acceleration to make dot phi smallest 
        dotphi = -sin(theta)*dtheta - k*cos(theta)*dtheta^2 - k*sin(theta);
        try 
            assert(dotphi < 0); % minimum dotphi should be negative
        catch
            disp(x)
            disp(dotphi);
        end
        u_s = 1;
    else
        u_s = u_r;
    end

    % simulate the control 
    x = next_state(x, u_s, dt);
end

% plot the safety index evolution 
plot(phis, 'LineWidth',3);
ylabel('safety index');
xlabel('time steps');


function x = next_state(x, u_s, dt)
    theta = x(1);
    dtheta = x(2);
    theta = clip(theta + dtheta*dt, 0, pi);
    dtheta = clip(dtheta + u_s*dt, -1, 1);
    
    x = [theta; dtheta];
end

function y = clip(x,bl,bu)
  % return bounded value clipped between bl and bu
  y=min(max(x,bl),bu);
end

function phi_val = phi(x)
    theta = x(1);
    dtheta = x(2);
%     k = 3.741078e-06;
%     k = 0.5324;
%     k = 10;
%     k = 3.3634;
%     k = 2.63;
    k = 2
    phi_val = cos(theta) - 0.5 - k * sin(theta) * dtheta;
end