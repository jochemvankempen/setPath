function noise = generateNoise(signalLength,fs,nTrials,noiseAmplitude,CATEGORY)
% genarates noise. Make sure to add the noise in frequency plane,
% especially the white noise. This is the only way 
switch CATEGORY
    case 'white'
        a=0;
        b=2*pi;
        noise = (a + (b-a).*rand(nTrials,signalLength)) ;
        
        noise = exp(1i*noise)*noiseAmplitude;
%         
%         
%         f = fft(noise')/size(noise,2);
%         f(1,:)=0;
%         noise = (ifft(f)*signalLength)';
    case 'gauss'
        noise = randn(nTrials,signalLength)*sqrt(0.1); 
        noise = exp(1i*noise)*noiseAmplitude;

end


if 0
    figure(1000),clf
    f = fft(noise(1,:))/size(noise(1,:),2);
%     meanf = mean(f,2);
    frequencies = linspace(0,fs/2,floor(length(f(1,:))/2));
    subplot(1,2,1)
    bar(frequencies,abs(f(1:floor(signalLength/2))).^2)

%     plot(abs(meanf))
%     ylim([-0.0002 2])
    xlim([-100 fs/2])
    axis square
    title(['barplot of ' CATEGORY ' noise'])
    subplot(1,2,2)
    hist(abs(f))
    axis square
    title(['hist of ' CATEGORY ' noise'])

end