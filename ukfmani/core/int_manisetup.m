% builds structures for fast pack and unpack, also for nested manifolds
%
% pack/unpack group (cell NxC) <-> (matrix NxG)
% pack/unpack alg   (cell NxC) <-> (matrix NxA)
%
% Emanuele Ruffaldi 2017 @ SSSA
function s = int_manisetup(group,alg,m)
if isfield(m,'models')
    s = [];
    for J=1:length(m.models)
        q = int_manisetup(group,alg,m.models{J});
        if isempty(s)
            s = q;
        else            
            s = [s, q];
        end
        group = group + m.models{J}.group;
        alg = alg + m.models{J}.alg;        
    end
else
    s =  struct('model',m,'group',group,'alg',alg);
    s.group = [group,group+m.group-1];
    s.alg = [alg,alg+m.alg-1];
end