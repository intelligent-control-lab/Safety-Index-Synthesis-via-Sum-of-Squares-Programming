%% this is design for the general safety index 
% find the sos polynomial: x^c1 + c1 x^(c1-1) + 1

% sos decomposition:
% Z = [1; x^{c1/2 - 1}; x^{c1/2}]

% decision variable is c1
options = optimoptions('fmincon','Display','iter','Algorithm','SQP','ConstraintTolerance',5e-15);
% options = optimoptions('fmincon','Display','iter','Algorithm','SQP');
obj = @(x)1;
disp(obj(100000));

xref = [9];

% just ensure c1 is positive, dont be too large
LB = [0];
UB = [10];

% find it 
A = [];
b = [];

% solve for solution 
[x, fval, exitflag,output] = fmincon(obj,xref,[],[],A,b,LB,UB,@(x)nonlcon_sdp(x),options);

fprintf("the computed solution is %d \n", floor(x));
% fprintf("the computed solution is %d \n", x);

function [c, ceq] = nonlcon_sdp(c1)
% first get the integer of c1 
% c1 = floor(c1);

% construct the Q matrix
if c1 == 2
    Q_sum = [1 c1/2;
             c1/2 1];
else
    Q_sum = [1 0 0;
     0 0 c1/2;
     0 c1/2 1];
end

% get the nonlinear constraints for SDP
c = [];
for i = 1:size(Q_sum,2)
    c = [c; -det(Q_sum(1:i,1:i))];
end
c = max(c);

% equality constraint
ceq = mod(c1,2);
end

