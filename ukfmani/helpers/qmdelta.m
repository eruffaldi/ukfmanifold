function w = qmdelta(X,Y)

w = getomega(qmult(X,qconj(Y))); % log(x*inv(y))

function w = getomega(q)

[v,phi] = qdecomp(q);
assert(isnan(phi) == 0)
w = v*phi;
w = w(:)';