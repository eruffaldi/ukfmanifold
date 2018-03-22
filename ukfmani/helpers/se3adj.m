% Adjoint of group (R,t) or of R as 4x4
%
% Convention: omega linear
function M = se3adj(R,t)

if nargin == 1
    if length(R) == 4
        t = R(1:3,4);
        R = R(1:3,1:3);
    else
        t = [0,0,0];
    end
end
M = blkdiag(R,R);
M(4:6,1:3) = skew(t)*R;
