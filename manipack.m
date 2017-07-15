% takes a model m and cell array and builds a paclged group
function y = manipack(m,c)

assert(isfield(m,'s'),'missing setup, use manisetup(m)');
y = zeros(size(c,1),m.group);
s = m.s;

for I=1:length(s)
    for J=1:size(c,1)
        y(J,s(I).group(1):s(I).group(2)) = s(I).model.pack(c{J,I});       
    end
end