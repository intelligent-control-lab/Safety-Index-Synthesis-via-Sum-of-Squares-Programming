%% the Unicycle Safety Index Design  
% ======= refute problem =======
% x - 5 > 0
% -x^2 + x > 0


%% the general safety index design 
clc
clear
% parameter settings
max_iter = 40;

% decision variable is c1
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',5e-13);
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',1e-15);
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','ConstraintTolerance',1e-15);

obj = @(x)1;
 
% decision variables [p1 p2]
% constraint the lower and upper bound
% LB = [zeros(1,4), -inf]; 
LB = [zeros(1,2)];
UB = inf*ones(1,2);

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
for k = 1
    for i = 1:max_iter
        % x = [p1, p2, p3, p4, p5]
        rng(seed);
        xref = [1*rand(1,2)];
%          xref  = [1*rand(1,4), 1*(-1 + 2*rand(1))];
        
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
% values = ones(1,size(xs,2));
% plot(xs, values, 'o', 'LineWidth', 0.5)
% xlim([0,60])


%% auxiliary functions
function [c, ceq] = nonlcon_sdp(k,x)
    p1 = x(1);
    p2 = x(2);

    constant = 1 - 5*p1;
    x_coe = p1+p2;
    x2_coe = -p2;

    % global variables = [1 x]
    % we first deal with upper right half
    Q = zeros(2,2);

    Q(1,1) = -constant;
    Q(2,2) = -x2_coe;
    Q(1,2) = -x_coe/2;

    % flip to the full matrix
    % make the matrix as symmetric
    % get the symmetric matrix Q 
    dim = 2;
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



function maxeig = max_eig(k,x)
    p1 = x(1);
    p2 = x(2);

    constant = 1 - 5*p1;
    x_coe = p1+p2;
    x2_coe = -p2;

    % global variables = [1 x]
    % we first deal with upper right half
    Q = zeros(2,2);

    Q(1,1) = -constant;
    Q(2,2) = -x2_coe;
    Q(1,2) = -x_coe/2;

    % flip to the full matrix
    % make the matrix as symmetric
    % get the symmetric matrix Q 
    dim = 2;
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