% SO3 as matrix
%
% Emanuele Ruffaldi 2017 @ SSSA
function m = makeRot()


m = [];
m.type = {'Rot'};
m.inv = @(X) flatten(unflatten(X)');
m.prod = @(X,Y) flatten(unflatten(X)*unflatten(Y));
m.log = @(X) so3log(unflatten(X));
m.exp = @(x) flatten(so3exp(x));
m.delta = @(X,Y) so3log(unflatten(X)*unflatten(Y)'); 
m.step = @(X,y) flatten(so3exp(y)*unflatten(X));
m.meancov = @manimeancov;
m.count = 1;
m.transport = @(X,t,Y) t;
m.group = 9; % as matrix
m.alg = 3;
m.pack = @(x) x(:)';
m.unpack = @(x) reshape(x,3,3);
m.islie = 1;
m.s = int_manisetup([],[],m);


function y = flatten(x)
y = x(:)';

function y = unflatten(x)
y = reshape(x,3,3);
