function [c, ceq] = nonlcon_sdp_eval(x)
    k = x(1);
    p1 = x(2);
    p2 = x(3);
    p3 = x(4);
    
    constant = 1 + p2 - p1*k + p3;
    x2_coe = -p3;
    y2_coe = -p2;
    x2y_coe = p1;

    % global variables = [1 x y xy]
    % we first deal with upper right half
    Q = zeros(3,3);

    Q(1,1) = -constant;
    Q(2,2) = -x2_coe;
    Q(3,3) = -y2_coe;
    Q(2,4) = -x2y_coe/2;

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