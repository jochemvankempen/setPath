function  [FreqFilt FreqPhasePower] = filtForMI(Freq2testLow,Freq2testHigh,data2filt,fs,powerPhase,normalize)

if Freq2testLow <0
    Freq2testLow=0;
end

if nargin < 6
    if Freq2testLow == 0
        Freq2testLow = Freq2testLow+1;%(1/(times2save(2)-times2save(1)));
    end
    normalize = 0;
elseif nargin == 6
    if Freq2testLow == 0
        Freq2testLow = Freq2testLow+1;%(1/(times2save(2)-times2save(1)));
    end
    
    disp('normalizing band passed power')
elseif nargin == 7
    if Freq2testLow == 0
        Freq2testLow = Freq2testLow+1;%(1/(times2save(2)-times2save(1)));
    end
    
%     Freq2testLow = Freq2testLow+(1/(times2save(2)-times2save(1)));
end


FreqFilt  = eegfilt(double(data2filt), fs, Freq2testLow, Freq2testHigh);
% FreqFilt  = eegfilt(double(FreqFilt), fs, [], Freq2testHigh);
% FreqFilt  = eegfilt(double(data2filt), fs, Freq2testLow, []);

if normalize == 1
    % normalize for higher frequency power (canolty 2006). subtract the
    % temporal mean and divide by the temporal power
%     sum(FreqFilt)/length(FreqFilt)
    
    FreqFilt = (FreqFilt - (sum(FreqFilt)/length(FreqFilt)))/std(FreqFilt) ;
end

if strcmpi(powerPhase, 'phase')
    FreqPhasePower = angle(hilbert(FreqFilt));
elseif strcmpi(powerPhase, 'power')
    FreqPhasePower = abs(hilbert(FreqFilt)).^2;
end


%%
