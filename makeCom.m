% combination
function m = makeCom(m1,m2)


g = m1.group;
a = m1.alg;

m = [];
m.inv = @(x) [m1.inv(x(1:g)), m2.inv(x(g+1:end))];
m.prod = @(x1,x2) [m1.prod(x1(1:g),x2(1:g)), m2.prod(x1(g+1:end),x2(g+1:end))];
if isfield(m1,'log') && isfield(m2,'log')
    m.log = @(x) [m1.log(x(1:g)), m2.log(x(g+1:end))];
    m.exp = @(x) [m1.exp(x(1:a)), m2.exp(x(a+1:end))];
end
m.delta = @(x,y) [m1.delta(x(1:g),y(1:g)), m2.delta(x(g+1:end),y(g+1:end))];
m.step =@(X,y) [m1.step(X(1:g),y(1:a)), m2.step(X(g+1:end),y(a+1:end))]; 
m.group = m1.group+m2.group;
m.alg = m1.alg+m2.alg;
m.count = m1.count+m2.count;
m.models = {m1,m2};
m.pack = []; %@(x) [m1.flatten(x{1}), m2.flatten(x{2})]; % this works with cells
m.unpack = []; %@(x) {m1.unflatten(x(1:g)), m1.unflatten(x(g+1:end))};

m.meancov = @meancovRn;





