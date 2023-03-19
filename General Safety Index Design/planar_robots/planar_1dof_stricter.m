%% the general safety index design 
clc
clear
% parameter settings
max_iter = 50;

% decision variable is c1
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',5e-13);
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',1e-15);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','ConstraintTolerance',1e-15);
obj = @(x)-x(1);
% obj = @(x)x(1);
% obj = @(x)1; 

% constraint the lower and upper bound
LB = [zeros(1,5) -inf];
% UB = inf*ones(1,6);
UB = [20 inf*ones(1,6)]; 

%% before the solution 
seed = 1;
% seed = 9;
xs = [];
full_xs = [];
seeds = [];
numericals = [];
exitflags = [];

%% solve for solution
for i = 1:max_iter
    % x = [k, p1, p2, p3, p4, q1]
    rng(seed);
%     xref  = [1000*rand(1), 1000*rand(1,4), 1000*(-1 + 2*rand(1,2))];
%     xref  = [10*rand(1), 1000*rand(1,4), 1000*(-1 + 2*rand(1,2))];
%     xref  = [100*rand(1), 100*rand(1,4), 100*(-1 + 2*rand(1))];
    xref  = [10, 100*rand(1,4), 100*(-1 + 2*rand(1))];
    
    % find it 
    A = [];
    b = [];
    
    % solve for solution 
    [x, fval, exitflag,output] = fmincon(obj,xref,[],[],A,b,LB,UB,@(x)nonlcon_sdp(x),options);

    % collect the valid solutions
    if exitflag >= 0
        % log the solution 
        % valid solution found 
        c = nonlcon_sdp(x);
%         numerical = c - 1e-12; % minus the 1e-12 offset from the constraints 
        numerical = c;
        if numerical < 1e-5 % sufficiently small numerical error
%         if 1
            % good solution 
            xs = [xs x(1)];
            full_xs = [full_xs; x];
            seeds = [seeds seed];
            numericals = [numericals numerical];
            exitflags = [exitflags exitflag];
        end
    end
%     break
    seed = seed + 1;
end

k = max(xs);
fprintf("the computed solution is %d \n", k);


function [c, ceq] = nonlcon_sdp(x)
    k = x(1);
    p1 = x(2);
    p2 = x(3);
    p3 = x(4);
    p4 = x(5);
    q1 = x(6);
    
    constant = 1 + p3 + p4 - q1 - sqrt(3)/2*p2;
    x_coe = p2-p3-p1*k;
    x2_coe = q1;
    y2_coe = q1;
    z2_coe = -p4;
    xz_coe = -p1;
    yz2_coe = -p1*k;

    % global variables = [1 x y z yz]
    % we first deal with upper right half
    Q = zeros(5,5);

    Q(1,1) = -constant;
    Q(1,2) = -x_coe/2;
    Q(2,2) = -x2_coe;
    Q(3,3) = -y2_coe;
    Q(4,4) = -z2_coe;
    Q(2,4) = -xz_coe/2;
    Q(4,5) = -yz2_coe/2;

    % flip to the full matrix
    % make the matrix as symmetric
    % get the symmetric matrix Q 
    dim = 5;
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
    c = [];
    for i = 1:size(Q,2)
        c = [c; -det(Q(1:i,1:i))];
    end
%     c = c + 1e-12;
    c = max(c);

    % equality constraint
    ceq = [];
end