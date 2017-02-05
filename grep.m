function [g, idx]  = grep(C, re)
% [g, idx]  = grep(C, re) similar to UNIX grep function
%   
% INPUTS: 
% C: a cell array of strings
% re: a matlab regular expression
% OUTPUTS:  
% g: the cell array of the strings matched by the regexp
% idx: the indces of the elements in C which matched
  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
%1

  
  
  [s, e] = regexpi(C, re);
  idx = ~cellfun('isempty', s);
  
  for i = 1:length(C)
    if isempty(s{i})
      g{i} = '';
    else
      g{i} = C{i}(s{i}(1):e{i}(1));
    end
  end
  
  