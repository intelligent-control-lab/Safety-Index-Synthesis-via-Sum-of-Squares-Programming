%% the verification that the SOS should be correct 
% sincos - k > 0
% sinx - 0.5 > 0
% 1 - sinx > 0
% sin^2 + cos^2 - 1 = 0 


%% the general safety index design 
clc
clear
% parameter settings
max_iter = 20;

% decision variable is c1
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','ConstraintTolerance',0);

obj = @(x)x(1); % minimize k 
 
% decision variables [k p1 p2 p3]
% constraint the lower and upper bound
LB = zeros(1,3); 
UB = inf*ones(1,3);
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
for k = 0.2:0.5:8
    for i = 1:max_iter
        % x = [k, p1, p2, p3]
        rng(seed);
        xref  = [100*rand(1,3)];
        
        % find it 
        A = [];
        b = [];
        
        % solve for solution 
        [x, fval, exitflag,output] = fmincon(obj,xref,[],[],A,b,LB,UB,@(x)nonlcon_sdp(k, x),options);
    
        % collect the valid solutions
        if exitflag >= 0
            % log the solution 
            % valid solution found 
            c = nonlcon_sdp(k, x);
    %         numerical = c - 1e-12; % minus the 1e-12 offset from the constraints 
            numerical = c;
            if numerical < 1e-15 % sufficiently small numerical error
    %         if 1
                % good solution 
                xs = [xs k];
                full_xs = [full_xs; x];
                seeds = [seeds seed];
                numericals = [numericals numerical];
                exitflags = [exitflags exitflag];
                break
            end
        end
    %     break
        seed = seed + 1;
    end
end


function [c, ceq] = nonlcon_sdp(k, x)
    p1 = x(1);
    p2 = x(2);
    p3 = x(3);
    
    constant = 1 - k*p1 + p3 + 1/3*p2;
    y_coe = -p1;
    x2_coe = -p2;
    y2_coe = -p3;
    xy2_coe = -k*p1;

    % global variables = [1 x y xy]
    % we first deal with upper right half
    Q = zeros(4,4);

    Q(1,1) = -constant;
    Q(1,3) = -y_coe/2;
    Q(2,2) = -x2_coe;
    Q(3,3) = -y2_coe;
    Q(3,4) = -xy2_coe/2;

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