function plotColorbar(varargin)
%
% INPUT:
% CATEGORY
%
% default is standard
%
% OUTPUT:
% a nice colobar
%
% Jochem van Kempen, 2014-apr-28

if nargin<1
    CATEGORY='standard';
else
    CATEGORY=varargin;
end

switch CATEGORY
    case 'standard'
        colormap default
        colorbar
    case 'circular'
        % for phase angles
%     case '
end