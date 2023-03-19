%% differential drive collision avoidance 

% according to hand derivation: 
% p0 = -F
% F(c1, c2)

% hyper parameter 
a_min = -1;
w_max = 1;
v_max = 1;
d_min = 1;

% decision variable is c1
options = optimoptions('fmincon','Display','iter','Algorithm','SQP','ConstraintTolerance',5e-13);
obj = @(x)1;

% x = [c2, p1, p2, p3, p4, p5, p6, t0, t1]
xref = zeros(1,9);
xref  = [10, zeros(1,8)];

% just ensure c1 is positive, dont be too large
LB = [zeros(1,7) -inf -inf];
UB = inf*ones(1,9);

% find it 
A = [];
b = [];

% solve for solution 
[x, fval, exitflag,output] = fmincon(obj,xref,[],[],A,b,LB,UB,@(x)nonlcon_sdp(x, a_min, w_max, v_max, d_min),options);

fprintf("the computed solution is %d \n", x(1));


function [c, ceq] = nonlcon_sdp(x,amin, wmax, vmax, dmin)
    % x = [c2, p1, p2, p3, p4, p5, p6, t0, t1]
    % get the variables
    c2 = x(1);
    p1 = x(2);
    p2 = x(3);
    p3 = x(4);
    p4 = x(5);
    p5 = x(6);
    p6 = x(7);
    t0 = x(8);
    t1 = x(9);

    % vector = [1 d v ca sa]
    % Q matrix is 5 by 5
    Q = zeros(5,5);
    Q(1,1) = -t1 + 1 + t0*dmin + p4*vmax; % constant 
    Q(4,4) = t1; % sin^2
    Q(5,5) = t1; % cos^2
    Q(1,5) = p6/2; % sin
    Q(1,3) = (p3 - p4)/2; % v
    Q(1,2) = (p2- t0)/2; % d
    Q(3,5) = -p1*c2*wmax / 2; % v * sin
    Q(1,4) = (c2*p1*amin + p5) / 2; % cos
    Q(3,4) = (p1 + t0*c2) / 2; % v * cos

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