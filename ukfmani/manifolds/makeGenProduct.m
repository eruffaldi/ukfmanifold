function mout = makeGenProduct(name,varargin)
% produce a new output file called: man##name.m
%
% that implements makeProduct of the above spec
%
% Supported:
% - Rot Rot1 SE3Mat SE3Quat Quat Quat1 Rn
% - Product ...
%
if nargin == 1
    mout = varargin{1};
    return;
end

q = makeProduct(varargin{:}); % then continue by using m.models m.types
ofile = ['man',name,'.m'];
if exist(ofile,'file')
    mout = eval([ofile(1:end-2),'()']);
    return;
end
q.output = ofile;
q.withlog = all(cellfun(@(x) isfield(x,'log'),varargin));
q = removefunction(q);
jq = jsonencode(q);
[pathstr,NAME,EXT] = fileparts(mfilename('fullpath'));
tmpfile = tempname();
disp(tmpfile)
fp = fopen(tmpfile,'w');
fprintf(fp,'%s',jq);
fclose(fp);
% store jq to tmpfile
system(['python ',pathstr,filesep,'makeGenProduct.py',' ',tmpfile]);
if exist(ofile,'file') == 0
    ofile = '';
    mout = [];
else
    mout = eval([ofile(1:end-2),'()']);
end


function x = removefunction(x)

if isstruct(x)
    f = fieldnames(x);
    for I=1:length(f)
        x.(f{I}) = removefunction(x.(f{I}));
    end
elseif iscell(x)
    for I=1:numel(x)
        x{I} = removefunction(x{I});
    end
elseif isa(x, 'function_handle')
    x = [];
end