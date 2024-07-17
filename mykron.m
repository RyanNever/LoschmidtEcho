function out = mykron(ma,mb,varargin)
% File: mykron.m
% Date: 17-Aug-98
% Author: I. Chuang <ike@isl.stanford.edu>
% kronecker product function which accepts multiple arguments.
if (length(varargin) == 0)
out = kron(ma,mb);
return;
else
out = mykron(kron(ma,mb),varargin{:});
return;
end
end

