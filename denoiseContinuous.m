function denoised = denoiseContinuous(x,fs, freqrange)

% denoisedX = denoiseContinuous(X,fs)
%
%Removes 60 Hz line noise and harmonics from a conintuous signal by dividing the signal into a set
%of overlapping windows and searching for peaks near multiples of 60Hz within each windowed segment.
% x - signal
%fs - sampling frequency
% 
% This function calls denoise60Hz

%C. Kovach 2008

x = double(x(:));

if nargin < 3
    freqrange = [59.5 60.5];
end

blocksize = round(300*fs./mean(freqrange)); %size of block in sampling points. It should be close to a multiple of line noise fundamental period.

%edge = 2000; % number of points on the edges to ignore due to edge effects
edge = 1000; % number of points on the edges to ignore due to edge effects
%avgWin = 2000; %Number of points to average with succeding trial to avoid discontinuities.
avgWin = ceil((blocksize - 2*edge - 1)./2);  %Number of points to average with succeding trial to avoid discontinuities.

N = length(x);

stepsize = blocksize - avgWin - 2*edge; 

denoised = zeros(size(x));

dnblock = denoise60Hz(x(1:blocksize,:),fs);

denoised(1:blocksize - edge,:) = dnblock(1:blocksize - edge,:);

steps = [stepsize:stepsize:N-blocksize];

%Using a function at least C2 continuous
avgFun = .5-.5*cos([0:avgWin-1]'./avgWin.*pi);

avgFun = avgFun(:,ones(size(x,2),1));

for s = steps
    
    dnblock = denoise60Hz(x(s:s+blocksize,:),fs,freqrange);
    
    %averaging together the ends
    denoised([s:s+avgWin-1]+edge,:) = (1-avgFun).*denoised([s:s+avgWin-1]+edge,:) + avgFun.*dnblock([1:avgWin]+edge,:);

    denoised([s:s + stepsize-1] + edge + avgWin,:) = dnblock([1:stepsize]+edge + avgWin,:);
end

dnblock = denoise60Hz(x(s+stepsize:end,:),fs);
denoised([s:s+avgWin-1]+edge + stepsize,:) = (1-avgFun).*denoised([s:s+avgWin-1]+edge + stepsize,:) + avgFun.*dnblock([1:avgWin]+edge,:);
denoised(s+edge+avgWin + stepsize:end,:) = dnblock(1+edge + avgWin:end,:);
    
