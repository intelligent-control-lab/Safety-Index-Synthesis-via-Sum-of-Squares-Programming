%% the Unicycle Safety Index Design 
% ======= original refute problem =======
% v - k - k*v*tan(alpha) > 0
% -v^2 + v > 0
% tan(alpha) > 0

% ======= refute problem after subsitution=======
% take v*tan(alpha) = y, v = x
% x - k - k*y > 0
% -x^2 + x > 0
% y > 0
% is empty


%% the general safety index design 
clc
clear
% parameter settings
max_iter = 10;

% decision variable is c1
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',5e-13);

obj = @(x)x(1);
 
% decision variables [k, p1 p2 p3 p4]
% constraint the lower and upper bound
% LB = [zeros(1,4), -inf]; 
LB = [zeros(1,5)];
UB = inf*ones(1,5);

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
% for k = 0.5:0.5:10
tic
for i = 1:max_iter
    % x = [k, p1, p2, p3, p4]
    rng(seed);
%     xref = [1*rand(1,5)];
    xref = [5, 1*rand(1,4)];
    
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
        maxeig = max_eig(x);
%             if numerical < 1e-10 % sufficiently small numerical error
        if numerical < 0 % strictly satisfy the constraint
            % good solution 
            xs = [xs x(1)];
            full_xs = [full_xs; x];
            max_eigs = [max_eigs maxeig];
            seeds = [seeds seed];
            numericals = [numericals numerical];
            exitflags = [exitflags exitflag];
        end
    end
    seed = seed + 1;
end
toc

fprintf('the minimum k for unicycle is designed as %d', min(xs));

%% auxiliary functions
function [c, ceq] = nonlcon_sdp(x)
    k = x(1);
    p1 = x(2);
    p2 = x(3);
    p3 = x(4);
    p4 = x(5);

    constant = 1 - k*p1;
    x_coe = p1+p2;
    y_coe = p3 - k*p4 - k*p1;
    x2_coe = -p2;
    y2_coe = -p4*k;
    xy_coe = p4;

    % global variables = [1 x y]
    % we first deal with upper right half
    Q = zeros(3,3);

    Q(1,1) = -constant;
    Q(2,2) = -x2_coe;
    Q(3,3) = -y2_coe;
    Q(1,2) = -x_coe/2;
    Q(1,3) = -y_coe/2;
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

    c = -eig(Q);
    c = max(c);

    % equality constraint
    ceq = [];
end



function maxeig = max_eig(x)
    k = x(1);
    p1 = x(2);
    p2 = x(3);
    p3 = x(4);
    p4 = x(5);

    constant = 1 - k*p1;
    x_coe = p1+p2;
    y_coe = p3 - k*p4 - k*p1;
    x2_coe = -p2;
    y2_coe = -p4*k;
    xy_coe = p4;

    % global variables = [1 x y]
    % we first deal with upper right half
    Q = zeros(3,3);

    Q(1,1) = -constant;
    Q(2,2) = -x2_coe;
    Q(3,3) = -y2_coe;
    Q(1,2) = -x_coe/2;
    Q(1,3) = -y_coe/2;
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

    % sdp problem 
    % get the nonlinear constraints for SDP
    c = -eig(Q);
    maxeig = max(c);
end