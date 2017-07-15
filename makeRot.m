% SO3 as matrix
function m = makeRot()


m = [];
m.inv = @(X) uflatten(unflat(X))';
m.prod = @(X,Y) uflatten(unflat(X)*unflat(Y));
m.log = @(X) so3log(unflatten(X));
m.exp = @(x) uflatten(so3exp(x));
m.delta = @(X,Y) so3log(unflatten(X)*unflatten(Y)'); 
m.step = @(X,y) uflatten(so3exp(y)*unflatten(X));
m.meancov = @manimeancov;
m.count = 1;
m.group = 9; % as matrix
m.alg = 3;
m.pack = @(x) x(:);
m.unpack = @(x) reshape(x,3,3);

function y = uflatten(x)
y = x(:);

function y = unflat(x)
y = reshape(x,3,3);
