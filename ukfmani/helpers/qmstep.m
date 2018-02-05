function q = mstep(X,y)

if norm(y) == 0
    q = X;
else
    q = qnorm(qmult(qomega2q(y(:)'),X));
end
