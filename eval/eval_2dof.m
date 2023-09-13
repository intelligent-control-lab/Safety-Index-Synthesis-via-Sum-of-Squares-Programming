clc
clear 

%% initialization  
max_iter = 10000000;
x = [pi/2,0,pi/2,0]; % initial state, theta = pi/2, dtheta = 0.
dt = 1e-6;
% k = 3.741078e-06;
% k = 0.5324;
% k = 10;
% k = 3.3634;
% k = 178;
k = 0.6;
save_freq = 20;

thetalim = [pi/3, 2*pi/3];
dthetalim = [-1,1];
ddthetalim = [-1, 1];


%% simulation 
phis = [];
for i = 1:max_iter
%     u_r = 2*rand(1) - 1; % random reference control
    u_r = [-1, -1]; % the maximum unsafe control 
    theta_1 = x(1);
    dtheta_1 = x(2);
    theta_2 = x(3);
    dtheta_2 = x(4);
    if mod(i, save_freq) == 0
        phis = [phis phi(x)];
    end
    if phi(x) > -1e-5 % due to discrete time system 
        % apply the maximum acceleration to make dot phi smallest 
        dotphi = -sin(theta_1)*dtheta_1 - k*cos(theta_1)*dtheta_1^2 - k*sin(theta_1) -sin(theta_2)*dtheta_2 - k*cos(theta_2)*dtheta_2^2 - k*sin(theta_2);
        try 
            assert(dotphi < 0); % minimum dotphi should be negative
        catch
            disp(x)
            disp(dotphi);
        end
        u_s = [1, 1];
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
    theta_1 = x(1);
    dtheta_1 = x(2);
    theta_2 = x(3);
    dtheta_2 = x(4);
    theta_1 = clip(theta_1 + dtheta_1*dt, pi/3, 2*pi/3);
    theta_2 = clip(theta_2 + dtheta_2*dt, pi/3, 2*pi/3);
    dtheta_1 = clip(dtheta_1 + u_s(1)*dt, -1, 1);
    dtheta_2 = clip(dtheta_2 + u_s(2)*dt, -1, 1);
    
    x = [theta_1, dtheta_1, theta_2, dtheta_2];
end

function y = clip(x,bl,bu)
  % return bounded value clipped between bl and bu
  y=min(max(x,bl),bu);
end

function phi_val = phi(x)
    theta_1 = x(1);
    dtheta_1 = x(2);
    theta_2 = x(3);
    dtheta_2 = x(4);
%     k = 3.741078e-06;
%     k = 0.5324;
%     k = 10;
%     k = 3.3634;
%     k = 2.63;
%     k = 1e-06;
%     k = 0.5;
    k = 0.6;
    phi_val = cos(theta_1) - k * sin(theta_1) * dtheta_1 + cos(theta_2) - k * sin(theta_2) * dtheta_2 - 1.0;
end