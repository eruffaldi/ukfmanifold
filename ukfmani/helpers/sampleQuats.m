function y = sampleQuats(n,muquat,Sigma)

C = cholcov(Sigma);
y = zeros(4,n);
Q = randn(3,n);
for i=1:n
    y(:,i) = mstep(muquat,C*Q(:,i));
end
y = y'; %n x 4
 

function q = mstep(X,y)

if norm(y) == 0
    q = X;
else
    q = qnorm(qmult(qomega2q(y(:)'),X));
end
