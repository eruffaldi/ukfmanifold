% takes a model m and cell array and builds a algebra
% Emanuele Ruffaldi 2017 @ SSSA
function y = manipackalg(m,c)

assert(isfield(m,'s'),'missing setup, use manisetup(m)');
y = zeros(size(c,1),m.alg);
s = m.s;

for I=1:length(s)
    for J=1:size(c,1)
        y(J,s(I).alg(1):s(I).alg(2)) = c{J,I};
    end
end