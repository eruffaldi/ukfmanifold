% we use this:
% x = {[1,2,3],[],1};
% y = cell(1,nargout(@max));
% [y{:}] = max(x{:})
%
%
% Also remember that anonymous can return multiple output%
% @(x,y) deal(x,y)
%
% Emanuele Ruffaldi 2017
function varargout = testanycall(varargin)

disp(sprintf('Nin: %d Nout:%d',nargin,nargout));

varargout = varargin;