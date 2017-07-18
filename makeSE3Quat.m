% SE3 to be tested, alternatively use makeCom(makeRot(),makeRn(3))
%
% Emanuele Ruffaldi 2017 @ SSSA
function m = makeSE3Quat()

m = [];
m.inv = @se3qinv;
m.prod = @se3qmul;
m.log = @se3qlog;
m.exp = @se3qexp;
m.delta = @se3qdelta;
m.step = @seqstep;
m.pack = @(x) x;
m.unpack = @(x) x;
m.group = 7; 
m.alg = 6;
m.count = 1;

function y = se3qdelta(x,y)
error('not implemented');

function y = se3qstep(x,y)
error('not implemented');

function y = se3qlog(x)
error('not implemented');

function y = se3qexp(x)
error('not implemented');

function y = se3qmul(x,y)

error('not implemented');

function y = se3qinv(x)

qc = qconj(x(1:4));
y = [qc, - qvqform(x(5:7))];

function u = ubuild2(quat,pos)

u = [quat(:);pos(:)]';

function q = getquat(x)

q = x(1:4);

function p = getpos(x)

p = x(5:7);

function u = getomega(x)

u = x(1:3);

function u = getvel(x)

u = x(4:6);