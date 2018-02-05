% replicates the model m0 n times
function m = makePower(m0,n)

m = makeProduct(repmat({m0},1,n));
m.type = {'Power',n,m0.type};