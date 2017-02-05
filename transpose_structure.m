function s = transpose_structure(s)
% S = TRANSPOSE_STRUCTURE(S) transposes all the data in the 
% structure S, e.g. 
% if s.x = ones(10,1), then 
% s = transpose_structure(s); % will end up with
% s.x = ones(1,10);
% Inputs:
%       s = Matlab structure
%
% Outputs:
%       s = structure with transposed fields
% 
% Usage: 
%       s = struct('x',ones(10,1));help transpose

%       ts = transpose_structure(s);
%       ts.x 
%
% Author: Alex Liberzon
% Copyright (c) 2013 alex.liberzon@gmail.com
% Turbulence Structure Laboratory, Tel Aviv University

f = fieldnames(s);
for i = 1:numel(s)
    for j = 1:numel(f)
        s(i).(f{j}) = s(i).(f{j}).'; 
    end
end
