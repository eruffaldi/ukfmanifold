% takes a model(C,G,A) m and a packed group x (matrix NxG) and builds a cell array
% (NxC)
function c = maniunpack(m,x)


assert(isfield(m,'s'),'missing setup, use manisetup(m)');
c = cell(size(x,1),m.count);
s = m.s;

for I=1:length(s)
    for J=1:size(x,1)
        c{J,I} = s(I).model.unpack(x(J,s(I).group(1):s(I).group(2)));
    end
end


