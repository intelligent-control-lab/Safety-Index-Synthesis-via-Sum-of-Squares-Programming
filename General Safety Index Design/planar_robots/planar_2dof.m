%% the general safety index design 


% decision variable is c1
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',1e-16);
obj = @(x)x(1);
% obj = @(x)1; 

% x = [k, p1, p2, p3, p4, q1, q2]
xref  = [10, 1000*rand(1,4), 1000*(-1 + 2*rand(1,2))];

% just ensure c1 is positive, dont be too large
LB = [zeros(1,5) -inf -inf];
UB = inf*ones(1,7);

% find it 
A = [];
b = [];

% solve for solution 
[x, fval, exitflag,output] = fmincon(obj,xref,[],[],A,b,LB,UB,@(x)nonlcon_sdp(x),options);

fprintf("the computed solution is %d \n", x(1));


function [c, ceq] = nonlcon_sdp(x)
    k = x(1);
    p1 = x(2);
    p2 = x(3);
    p3 = x(4);
    p4 = x(5);
    q1 = x(6);
    q2 = x(7);
    
    constant = 1 - p1*k + p3 + p4 - q1 - 0.5*q2;
    x_coe = p2-p3;
    y_coe = q2;
    x2_coe = q1;
    y2_coe = q1;
    z2_coe = -p4;
    xz_coe = -(k*q2 + p1);
    yz2_coe = -p1*k;

    % global variables = [1 x y z yz]
    % we first deal with upper right half
    Q = zeros(5,5);
    Q(1,1) = constant;
    Q(1,2) = x_coe/2;
    Q(1,3) = y_coe/2;
    Q(2,2) = x2_coe;
    Q(3,3) = y2_coe;
    Q(4,4) = z2_coe;
    Q(2,4) = xz_coe/2;
    Q(4,5) = yz2_coe/2;

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
    Q = -Q;

    % sdp problem 
    % get the nonlinear constraints for SDP
    c = [];
    for i = 1:size(Q,2)
        c = [c; -det(Q(1:i,1:i))];
    end
    c = max(c);

    % equality constraint
    ceq = [];
end