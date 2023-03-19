%% the general safety index design 
clc
clear
% parameter settings
max_iter = 30;

% decision variable is c1
options = optimoptions('fmincon','Display','off','Algorithm','interior-point','ConstraintTolerance',5e-13);
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',1e-15);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','ConstraintTolerance',1e-15);
% obj = @(x)-x(1);
obj = @(x)x(1);
% obj = @(x)1; 

% constraint the lower and upper bound
LB = [zeros(1,9) -inf, -inf]; %[p1 p2 p3 p4 q1]
UB = inf*ones(1,11); %[p1 p2 p3 p4 q1]
% UB = [20 inf*ones(1,6)]; 

%% before the solution 
seed = 1;
% seed = 9;
xs = [];
full_xs = [];
seeds = [];
numericals = [];
exitflags = [];
max_eigs = [];

%% solve for solution
for k = 2:0.5:5
    for i = 1:max_iter
        % x = [p1, p2, p3, p4, q1]
        rng(seed);
        xref  = [100*rand(1,9), 100*(-1 + 2*rand(1)), 100*(-1 + 2*rand(1))];
        
        % find it 
        A = [];
        b = [];
        
        % solve for solution 
        [x, fval, exitflag,output] = fmincon(obj,xref,[],[],A,b,LB,UB,@(x)nonlcon_sdp(k, x),options);
        disp("k=")
        disp(k)
        disp("exitflag=")
        disp(exitflag)
        disp("x=")
        disp(x)
        % collect the valid solutions
        if exitflag >= 0
            % log the solution 
            % valid solution found 
            c = nonlcon_sdp(k, x);
    %         numerical = c - 1e-12; % minus the 1e-12 offset from the constraints 
            numerical = c;
            
            if numerical < 0 
    %         if 1
                % good solution 
                xs = [xs k];
                full_xs = [full_xs; x];
                seeds = [seeds seed];
                numericals = [numericals numerical];
                exitflags = [exitflags exitflag];
                break % no need to continue 
            end
        end
    %     break
        seed = seed + 1;
    end
end

k = min(xs);
fprintf("the computed solution is %d \n", k);


function [c, ceq] = nonlcon_sdp(k, x)
%     k = x(1);
%     p1 = x(2);
%     p2 = x(3);
%     p3 = x(4);
%     p4 = x(5);
%     q1 = x(6);

    p1 = x(1);
    %link1
    p2 = x(2);
    p3 = x(3);
    p4 = x(4);
    p5 = x(5);
    %link2
    p6 = x(6);
    p7 = x(7);
    p8 = x(8);
    p9 = x(9);
    %link1 equality
    q1 = x(10);
    %link2 equality
    q2 = x(11);
    
    constant = 1 + p3 + p4 - q1 - sqrt(3)/2*p2 + p7 + p8 - q2 - sqrt(2.5)/2*p6;
    % link1
    a_coe_1 = p2-p3-p1*k;
    a2_coe_1 = q1; % sin^2
    b2_coe_1 = q1; % cos^2
    y2_coe_1 = -p4;
    z_coe_1 = p5; % dtheta^2
    z2_coe_1 = -p5;
    ay_coe_1 = -p1; % sin * dtheta
    bz_coe_1 = -p1*k; % cos * dtheta^2
    % link2
    a_coe_2 = p6-p7-p1*k;
    a2_coe_2 = q2; % sin^2
    b2_coe_2 = q2; % cos^2
    y2_coe_2 = -p8;
    z_coe_2 = p9; % dtheta^2
    z2_coe_2 = -p9;
    ay_coe_2 = -p1; % sin * dtheta
    bz_coe_2 = -p1*k; % cos * dtheta^2

    % global variables = [1 y_1 z_1 a_1 b_1 y_2 z_2 a_2 b_2] 
    % y=theta
    % z=dtheta^2
    % a=sin(theta)
    % b=cos(theta)
    
    % we first deal with upper right half
    Q = zeros(9,9);

    Q(1,1) = -constant;
    % link1
    Q(1,4) = -a_coe_1/2;
    Q(4,4) = -a2_coe_1;
    Q(5,5) = -b2_coe_1;
    Q(3,3) = -z2_coe_1;
    Q(2,2) = -y2_coe_1;
    Q(1,3) = -z_coe_1/2;
    Q(2,4) = -ay_coe_1/2;
    Q(3,5) = -bz_coe_1/2;
    % link2
    Q(1,8) = -a_coe_2/2;
    Q(8,8) = -a2_coe_2;
    Q(9,9) = -b2_coe_2;
    Q(7,7) = -z2_coe_2;
    Q(6,6) = -y2_coe_2;
    Q(1,7) = -z_coe_2/2;
    Q(6,8) = -ay_coe_2/2;
    Q(7,9) = -bz_coe_2/2;

    % flip to the full matrix
    % make the matrix as symmetric
    % get the symmetric matrix Q 
    dim = 9;
    for i = 1:dim
        for j = 1:i-1
            assert(j < i);
            Q(i,j) = Q(j,i);
        end
    end

    % take the negative
%     Q = -Q;

    % sdp problem 
    % get the nonlinear constraints for SDP
%     c = [];
%     for i = 1:size(Q,2)
%         c = [c; -det(Q(1:i,1:i))];
%     end
%     c = max(c);
% 
%     % equality constraint
%     ceq = [];

    c = -eig(Q);
    c = max(c);

    % equality constraint
    ceq = [];
end