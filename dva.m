function [d,mn,n] = dva(vecSet,varargin)
%
% function [d,mn,n] = dva(vecSet)
%
% Compute angular dispersion as the norm of the sum of norm'd vectors.
%

sz = size(vecSet);
nVec = sz(2);
u = zeros(sz);
n = zeros(nVec,1);
for i=1:nVec
    n(i) = norm(vecSet(:,i));
    u(:,i) = vecSet(:,i)./n(i); % divide each vector by its norm
end

d = 1 - norm(sum(u,2)) ./ nVec;
mn = mean(n);
