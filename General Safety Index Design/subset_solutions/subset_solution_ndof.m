%% the general safety index design 
clc
clear
% how many links of robot arms
global num_links;
% num_links = 6;
num_links = 12;
% num_links = 24;
global limitation_links;
% limitation_links = [sqrt(3)/2, sqrt(2.5)/2, sqrt(2.7)/2, sqrt(2.3)/2, sqrt(2.8)/2, sqrt(2.6)/2];
limitation_links = [sqrt(3)/2, sqrt(2.5)/2, sqrt(2.7)/2, sqrt(2.3)/2, sqrt(2.8)/2, sqrt(2.6)/2, sqrt(3)/2, sqrt(2.5)/2, sqrt(2.7)/2, sqrt(2.3)/2, sqrt(2.8)/2, sqrt(2.6)/2];
% limitation_links = [sqrt(3)/2, sqrt(2.5)/2, sqrt(2.7)/2, sqrt(2.3)/2, sqrt(2.8)/2, ...
%                     sqrt(2.6)/2, sqrt(3)/2, sqrt(2.5)/2, sqrt(2.7)/2, sqrt(2.3)/2, ...
%                     sqrt(2.8)/2, sqrt(2.6)/2, sqrt(3)/2, sqrt(2.5)/2, sqrt(2.7)/2, ...
%                     sqrt(2.3)/2, sqrt(2.8)/2, sqrt(2.6)/2, sqrt(3)/2, sqrt(2.5)/2, ...
%                     sqrt(2.7)/2, sqrt(2.3)/2, sqrt(2.8)/2, sqrt(2.6)/2];


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
LB = [zeros(1, 1 + 4 * num_links), -inf * ones(1, num_links)]; %[p1 p2 p3 p4 q1]
UB = inf*ones(1, 1 + 4 * num_links + num_links); %[p1 p2 p3 p4 q1]
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
% 3.5 is the solution for 6
for k = 4.5:0.5:10
% for k = 4.5
% for k = 4
    for i = 1:max_iter
        % x = [p1, p2, p3, p4, q1]
        rng(seed);
        xref  = [100*rand(1,1 + 4 * num_links), 100*(-1 + 2*rand(1, num_links))];
        
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

% k = max(xs);
% fprintf("the computed solution is %d \n", k);


function [c, ceq] = nonlcon_sdp(k, x)
    global num_links;
    global limitation_links;
    p1 = x(1);
    constant = 1;
    Q = zeros(1+4*num_links, 1+4*num_links);
        
    
    for i = 1:num_links
        %link-i inequality
        p2 = x((i-1)*4+2);
        p3 = x((i-1)*4+3);
        p4 = x((i-1)*4+4);
        p5 = x((i-1)*4+5);
        %link-i equality
        q1 = x(1+4*num_links+i);
        
        constant = constant + p3 + p4 - q1 - limitation_links(i)*p2;
        a_coe = p2-p3-p1*k;
        a2_coe = q1; % sin^2
        b2_coe = q1; % cos^2
        y2_coe = -p4;
        z_coe = p5; % dtheta^2
        z2_coe = -p5;
        ay_coe = -p1; % sin * dtheta
        bz_coe = -p1*k; % cos * dtheta^2
        
        % link-i           global variables = [1 ....  y_i z_i a_i b_i ....] 
        Q(1,4+(i-1)*4) = -a_coe/2;
        Q(4+(i-1)*4,4+(i-1)*4) = -a2_coe;
        Q(5+(i-1)*4,5+(i-1)*4) = -b2_coe;
        Q(3+(i-1)*4,3+(i-1)*4) = -z2_coe;
        Q(2+(i-1)*4,2+(i-1)*4) = -y2_coe;
        Q(1,3+(i-1)*4) = -z_coe/2;
        Q(2+(i-1)*4,4+(i-1)*4) = -ay_coe/2;
        Q(3+(i-1)*4,5+(i-1)*4) = -bz_coe/2;
    end
    
    Q(1,1) = -constant;
    % flip to the full matrix
    % make the matrix as symmetric
    % get the symmetric matrix Q 
    dim = 1+4*num_links;
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