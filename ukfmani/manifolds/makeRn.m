% euclidean
%
% Emanuele Ruffaldi 2017 @ SSSA
function m = makeRn(n)


m = [];
m.type = {'Rn',n};
m.fastinv = 1;
m.inv = @(x) -x;
%m.prod = @(x,y) x+y;
m.prod = @makecomp;
m.log = @(x) x;
m.exp = @(x) x(:)';
m.meancov = @meancov;
m.delta = @(x,y) x-y;
m.step = @(x,y) x+y;
m.group = n;
m.alg = n;
m.transport = @(X,t,Y) t;
m.count = 1;
%m.pack = @(x) x;
m.pack = @mpack;
m.unpack = @(x) x;
m.islie = 1;
m.s = int_manisetup([],[],m);


% added test functions mpack, makecomp
function p = mpack(x)
if iscell(x)
    p = reshape(x{1},[],1);
else
    p = x;
end

function c = makecomp(x, y)
c = x + y;

function [m,C] = meancov(x)

m = mean(x);
C = cov(x);
