% SE3 in 4x4 matrix
%
% Emanuele Ruffaldi 2017 @ SSSA
function m = makeSE3Mat()

m1= makeRot();
m = [];
m.inv = @se3inv;
m.prod = @(x,y) mflat(munflat(x)*munflat(y));
m.m1 = makeRot();
m.log = @(x) se3log(m1,x);
m.exp = @(x) se3exp(m1,x);
m.delta = @(x,y) se3delta(m1,x,y);
m.pack = @mflat;
m.unpack = @munflat;
m.step = @(X,y) mflat(munflat(m.exp(y))*munflat(X));
m.group = 16;
m.alg = 6;
m.count = 1;

function y = se3log(m1,x)
y = [m1.log(getrot(x)), getpos(x)]; % not exact

function y = se3exp(m1,x)

y = mflat(ubuild2(m1.exp(getomega(x)),x(4:6))); % not exact


function r = se3delta(m1,x,y)

r =  [m1.delta(getrot(x),getrot(y)), getpos(x)-getpos(y)];

function u = mflat(x)

u = x(:)';

function u = munflat(x)

u = reshape(x,4,4);

function y = se3inv(x)

R = getrot(x);
y = eye(4);
y(1:3,1:3) = R';
y(1:3,4) = -y(1:3,1:3)*getpos(x)';
y = mflat(y);

function u = ubuild2(rot,pos)

u = eye(4);
u(1:3,1:3) = reshape(rot,3,3);
u(1:3,4) = pos;

function R = getrot(x)

M = munflat(x);
R = M(1:3,1:3);

function p = getpos(x)

p = x(13:15);
p = p(:)';

function u = getomega(x)

u = x(1:3);

function u = getvel(x)

u = x(4:6);