%% the verification that the SOS should be correct 
% sincos - k > 0
% sinx - 0.5 > 0
% 1 - sinx > 0
% sin^2 + cos^2 - 1 = 0 


%% the general safety index design 
clc
clear
% parameter settings
max_iter = 50;

% decision variable is c1
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',5e-13);
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',1e-15);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','ConstraintTolerance',1e-15);

obj = @(x)x(1); % minimize k 
 
% decision variables [k p1 p2 p3 q1]
% constraint the lower and upper bound
LB = [zeros(1,4) -inf]; 
UB = inf*ones(1,5);
% UB = [20 inf*ones(1,6)]; 

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
    % x = [k, p1, p2, p3, q1]
    rng(seed);
    xref  = [100, 100*rand(1,3), 100*(-1 + 2*rand(1))];
    
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
        if numerical < 1e-6 % sufficiently small numerical error
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
    q1 = x(5);
    
    constant = 1 - k*p1 - 0.5*p2 + p3 - q1;
    x_coe = p2-p3;
    x2_coe = q1;
    y2_coe = q1;
    xy_coe = p1;

    % global variables = [1 x y]
    % we first deal with upper right half
    Q = zeros(3,3);

    Q(1,1) = -constant;
    Q(1,2) = -x_coe/2;
    Q(2,2) = -x2_coe;
    Q(3,3) = -y2_coe;
    Q(2,3) = -xy_coe/2;

    % flip to the full matrix
    % make the matrix as symmetric
    % get the symmetric matrix Q 
    dim = 3;
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