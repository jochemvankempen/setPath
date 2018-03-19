function [swv,swn] = fast_sw_dva(vecSet,winSize)
%
% function [swv,swn] = fast_sw_dva(vecSet,winSize)
%
% Computes dva on the vector set in a sliding window of size winSize.
% vecSet should be [dim x time]
% Best if winSize is an odd number.
% Fast (0/1) means use fast version (why not?).
%
% AS Apr 2014
%
% 14 May 2014
% Added runNorm, so that running mean norm is also computed incrementally.
%

winSize = winSize-1; % makes things simpler. You end up with winSize.

sz = size(vecSet);
nTP = sz(2);
halfWin = floor(winSize/2);
swv = zeros(1,sz(2));
swn = zeros(1,sz(2));
n = zeros(1,sz(2));

% compute the NORM
for j=1:nTP
    n(j) = norm(vecSet(:,j));
end

% divide each vector by its norm
vecSet = bsxfun(@rdivide,vecSet,n);
nVec = sum(~isnan(nansum(vecSet,1)));

% start the incremental computation
thisWin = [max(1-halfWin,1):min(1+halfWin,nTP)];
nInWin = zeros(1,nTP);
nInWin(1) = length(thisWin);
runSum = nansum(vecSet(:,thisWin),2);
swv(1) = 1 - norm(runSum) ./ nInWin(1);
swn(1) = sum(n(thisWin));
runNorm = swn(1);


for j=2:nTP
    thisWin = [max(j-halfWin,1):min(j+halfWin,nTP)];
    nInWin(j) = length(thisWin);
    if thisWin(1)>1
        runSum = runSum - vecSet(:,thisWin(1)-1);
        runNorm = runNorm - n(thisWin(1)-1);
    end
    if thisWin(nInWin(j))<nTP
        runSum = runSum + vecSet(:,thisWin(nInWin(j)));
        runNorm = runNorm + n(thisWin(nInWin(j)));
    end
    swv(j) = 1 - norm(runSum) ./ nInWin(j);
    swn(j) = runNorm; %mean(n(thisWin));
end
swn = swn ./ nInWin;
