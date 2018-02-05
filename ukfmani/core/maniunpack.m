% takes a model(C,G,A) m and a packed group x (matrix NxG) and builds a cell array
% (NxC)
% Emanuele Ruffaldi 2017 @ SSSA
function c = maniunpack(m,x)


assert(isfield(m,'s'),'missing setup, use manisetup(m)');
if length(m.s) == 1
    c = m.unpack(x')'; % transposed because of convention 
else
    c = cell(size(x,1),m.count);
    s = m.s;

    for I=1:length(s)
        sm = s(I).model.unpack;
        for J=1:size(x,1)
            c{J,I} = sm(x(J,s(I).group(1):s(I).group(2))');
        end
    end
end
