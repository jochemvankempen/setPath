function denoised = denoise60hz(x,fs,frange)

%60 Hz denoising based on identifying periodic signals.
%
% x - a matrix of data with rows representing time and columns representing trials
% fs - sampling frequency
%frange - range of frequencies to search for line noise. If not specified,
%       then it is [59.5 61.5].

if nargin < 3
    frange = [59.5 61.5]; %Range of frequencies over which to look
end
step = .01; %Size of frequency step within the range
N = size(x,1);




nfrange = frange./(fs./2); %Converting frequencies to normalized 
nstep = step./(fs./2);


testfreq = nfrange(1):nstep:nfrange(2); % The script looks at power in fundamental and harmonics of these frequencies

cwidths = testfreq.*N./2; %Frequency steps in terms of sampling units

%The following constructs a set of dirac combs which correspond to 
% periodic signals with fundamental frequencies in testfreq

Combs = zeros(size(x,1),length(cwidths));

for i = 1:length(cwidths)
    
    Combs(round([cwidths(i):cwidths(i):N./2]+1),i) = 1;
    Combs(round(end - [cwidths(i):cwidths(i):N./2] + 1),i) = 1;
        
end

pow = abs(fft(x)).^2'*Combs;  %Amount of power explained by each comb 

[mx,maxcomb] = max(pow'); %Locating maximum for each trial

noiseft = fft(x).*Combs(:,maxcomb); %Noise is extracted by multiplying the fft of the signal by the best comb
                                   
noise = ifft(noiseft);

denoised = x - noise; %denoised signal
