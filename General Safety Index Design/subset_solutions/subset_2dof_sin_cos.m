%% the verification that the SOS should be correct 
% refute problem: to find the upper bound of y + xy^2
% -y - kxy^2 - k > 0 
% 1/3 - x^2 > 0 
% 1 - y^2 > 0 
% ======= equivalency =======
% -y -kxz - k > 0
% 1/3 - x^2 > 0
% 1 - y^2 > 0
% -z^2 + z > 0


%% the general safety index design 
clc
clear
% parameter settings
max_iter = 5;

% decision variable is c1
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',5e-13);
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',1e-15);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','ConstraintTolerance',1e-15);

obj = @(x)1;
 
% decision variables [p1 p2 p3 p4]
% constraint the lower and upper bound
LB = [zeros(1,7)]; 
UB = inf*ones(1,7);

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
for k = 0:0.1:1
    for i = 1:max_iter
        % x = [p1, p2, p3, p4]
        rng(seed);
        xref = [10*rand(1,7)];
        
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
            maxeig = max_eig(k,x);
%             if numerical < 1e-10 % sufficiently small numerical error
            if numerical < 0 % strictly satisfy the constraint
                % good solution 
                xs = [xs k];
                full_xs = [full_xs; x];
                max_eigs = [max_eigs maxeig];
                seeds = [seeds seed];
                numericals = [numericals numerical];
                exitflags = [exitflags exitflag];
                break % found the good k, just skip to next k 
            end
        end
        seed = seed + 1;
    end
end

%% draw 
% now draw all the feasible k solutions: 
values = ones(1,size(xs,2));
plot(xs, values, 'o', 'LineWidth', 0.5)
xlim([0,10])
% k = max(xs);
% fprintf("the computed solution is %d \n", k);

%% auxiliary functions
function [c, ceq] = nonlcon_sdp(k,x)
    p1 = x(1);
    p2 = x(2);
    p3 = x(3);
    p4 = x(4);
    p5 = x(5);
    p6 = x(6);
    p7 = x(7);

    constant = 1 - k*p1 + 1/3*p2 + p3 - k*p5 + 1/3*p6;
    x2_coe_1 = -p2;
    y2_coe_1 = -p3;
    z2_coe_1 = -p4;
    y_coe_1 = -p1;
    z_coe_1 = p4;
    xz_coe_1 = -p1*k;
    
    x2_coe_2 = -p5;
    y2_coe_2 = -p6;
    z2_coe_2 = -p7;
    y_coe_2 = -p1;
    z_coe_2 = p7;
    xz_coe_2 = -p1*k;
    
    

    % global variables = [1 x_1 y_1 z_1 x_2 y_2 z_2]
    % we first deal with upper right half
    Q = zeros(7,7);

    Q(1,1) = -constant;
    
    Q(2,2) = -x2_coe_1;
    Q(3,3) = -y2_coe_1;
    Q(4,4) = -z2_coe_1;
    Q(1,3) = -y_coe_1/2;
    Q(1,4) = -z_coe_1/2;
    Q(2,4) = -xz_coe_1/2;
 
    Q(5,5) = -x2_coe_2;
    Q(6,6) = -y2_coe_2;
    Q(7,7) = -z2_coe_2;
    Q(1,6) = -y_coe_2/2;
    Q(1,7) = -z_coe_2/2;
    Q(5,7) = -xz_coe_2/2;

    % flip to the full matrix
    % make the matrix as symmetric
    % get the symmetric matrix Q 
    dim = 7;
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
    p5 = x(5);
    p6 = x(6);
    p7 = x(7);

    constant = 1 - k*p1 + 1/3*p2 + p3 - k*p5 + 1/3*p6;
    x2_coe_1 = -p2;
    y2_coe_1 = -p3;
    z2_coe_1 = -p4;
    y_coe_1 = -p1;
    z_coe_1 = p4;
    xz_coe_1 = -p1*k;
    
    x2_coe_2 = -p5;
    y2_coe_2 = -p6;
    z2_coe_2 = -p7;
    y_coe_2 = -p1;
    z_coe_2 = p7;
    xz_coe_2 = -p1*k;
    
    

    % global variables = [1 x_1 y_1 z_1 x_2 y_2 z_2]
    % we first deal with upper right half
    Q = zeros(7,7);

    Q(1,1) = -constant;
    
    Q(2,2) = -x2_coe_1;
    Q(3,3) = -y2_coe_1;
    Q(4,4) = -z2_coe_1;
    Q(1,3) = -y_coe_1/2;
    Q(1,4) = -z_coe_1/2;
    Q(2,4) = -xz_coe_1/2;
 
    Q(5,5) = -x2_coe_2;
    Q(6,6) = -y2_coe_2;
    Q(7,7) = -z2_coe_2;
    Q(1,6) = -y_coe_2/2;
    Q(1,7) = -z_coe_2/2;
    Q(5,7) = -xz_coe_2/2;

    % flip to the full matrix
    % make the matrix as symmetric
    % get the symmetric matrix Q 
    dim = 7;
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