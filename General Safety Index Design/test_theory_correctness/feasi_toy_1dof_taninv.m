%% the verification that the SOS should be correct 
% refute problem: to find the upper bound of y + xy^2
% y + xy^2 - k > 0 (Xy - k > 0)
% 1 - x^2 > 0 (-X^2 + X > 0)
% 1 - y^2 > 0 (1 - y^2)
% ======= equivalency =======
% y + xz - k > 0
% 1 - x^2 > 0
% 1 - y^2 > 0
% -z^2 + z > 0


%% the general safety index design 
clc
clear
% parameter settings
max_iter = 10;

% decision variable is c1
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',5e-13);
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',1e-15);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','ConstraintTolerance',1e-15);

obj = @(x)1;
 
% decision variables [p1 p2 p3 p4]
% constraint the lower and upper bound
LB = [zeros(1,4)]; 
UB = inf*ones(1,4);

%% before the solution 
seed = 1;
k = 1;
% seed = 9;
xs = [];
full_xs = [];
seeds = [];
numericals = [];
exitflags = [];
max_eigs = [];

%% solve for solution
for i = 1:max_iter
    % x = [p1, p2, p3, p4]
    rng(seed);
    xref = [10*rand(1,4)];
    
    % find it 
    A = [];
    b = [];
    
    % solve for solution 
    [x, fval, exitflag,output] = fmincon(obj,xref,[],[],A,b,LB,UB,@(x)nonlcon_sdp(k,x),options);

    % collect the valid solutions
    if exitflag >= 0
        % log the solution 
        % valid solution found 
        c = nonlcon_sdp(k,x);
%         numerical = c - 1e-12; % minus the 1e-12 offset from the constraints 
        numerical = c;
%         if numerical < 1e-10 % sufficiently small numerical error
        maxeig = max_eig(k,x);
        if 1
            % good solution 
            xs = [xs x(1)];
            full_xs = [full_xs; x];
            max_eigs = [max_eigs maxeig];
            seeds = [seeds seed];
            numericals = [numericals numerical];
            exitflags = [exitflags exitflag];
        end
    end
%     break
    seed = seed + 1;
end

% k = max(xs);
% fprintf("the computed solution is %d \n", k);


function [c, ceq] = nonlcon_sdp(k,x)
    p1 = x(1);
    p2 = x(2);
    p3 = x(3);
    p4 = x(4);

    constant = 1 - k*p1 + p2 + p3;
    x2_coe = -p2;
    y2_coe = -p3;
    z2_coe = -p4;
    y_coe = p1;
    z_coe = p4;
    xz_coe = p1;

    

    % global variables = [1 x y z]
    % we first deal with upper right half
    Q = zeros(4,4);

    Q(1,1) = -constant;
    Q(2,2) = -x2_coe;
    Q(3,3) = -y2_coe;
    Q(4,4) = -z2_coe;
    Q(1,3) = -y_coe/2;
    Q(1,4) = -z_coe/2;
    Q(2,4) = -xz_coe/2;

    % flip to the full matrix
    % make the matrix as symmetric
    % get the symmetric matrix Q 
    dim = 4;
    for i = 1:dim
        for j = 1:i-1
            assert(j < i);
            Q(i,j) = Q(j,i);
        end
    end

    % sdp problem 
    % get the nonlinear constraints for SDP
%     c = [];
%     for i = 1:size(Q,2)
%         c = [c; -det(Q(1:i,1:i))];
%     end
% %     c = c + 1e-12;
%     c = max(c);

    c = -eig(Q);
    c = max(c);

    % equality constraint
    ceq = [];
end



function maxeig = max_eig(k,x)
    p1 = x(1);
    p2 = x(2);
    p3 = x(3);
    p4 = x(4);

    constant = 1 - k*p1 + p2 + p3;
    x2_coe = -p2;
    y2_coe = -p3;
    z2_coe = -p4;
    y_coe = p1;
    z_coe = p4;
    xz_coe = p1;

    

    % global variables = [1 x y z]
    % we first deal with upper right half
    Q = zeros(4,4);

    Q(1,1) = -constant;
    Q(2,2) = -x2_coe;
    Q(3,3) = -y2_coe;
    Q(4,4) = -z2_coe;
    Q(1,3) = -y_coe/2;
    Q(1,4) = -z_coe/2;
    Q(2,4) = -xz_coe/2;

    % flip to the full matrix
    % make the matrix as symmetric
    % get the symmetric matrix Q 
    dim = 4;
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