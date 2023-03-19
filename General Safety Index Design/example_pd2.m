%% this is design for the general safety index 
% find the sos polynomial: x^c1 + c1 x^(c1-2) + 1

% sos decomposition:
% Z = [1; x^{c1/2 - 1}; x^{c1/2}]

% decision variable is c1
options = optimoptions('fmincon','Display','iter','Algorithm','SQP','ConstraintTolerance',5e-15);
% options = optimoptions('fmincon','Display','iter','Algorithm','SQP');
obj = @(x)1;
disp(obj(100000));

xref = [10.35];

% just ensure c1 is positive, dont be too large
LB = [0];
UB = [100];

% find it 
A = [];
b = [];

% solve for solution 
[x, fval, exitflag,output] = fmincon(obj,xref,[],[],A,b,LB,UB,@(x)nonlcon_sdp(x),options);

fprintf("the computed solution is %d \n", x);


function [c, ceq] = nonlcon_sdp(c1)
% construct the Q matrix

Q_sum = [1 0 0;
        0 c1 0;
        0 0 1];

% get the nonlinear constraints for SDP
c = [];
for i = 1:size(Q_sum,2)
    c = [c; -det(Q_sum(1:i,1:i))];
end
c = max(c);

% equality constraint
ceq = mod(c1,2);
end

