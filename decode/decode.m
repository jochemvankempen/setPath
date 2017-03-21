function decode(savefilename,data,trialLabel,decodeDir,times2downsample,freq,DEC,trialTypes)

if exist([decodeDir 'dTF_' savefilename '.mat' ],'file'), disp([savefilename]), return, end
disp(['working on: ' 'dTF_' savefilename '.mat'])

switch DEC.dimension
    case {'eachElectrodeEachTime','multiElectrodeEachTime'};
        decode_eachTime_parallel(data, savefilename, trialLabel,decodeDir,times2downsample,freq,DEC,trialTypes)
    case {'eachElectrodeEachFreq','multiElectrodeEachFreq'};
        decode_eachFreq_parallel(data, savefilename, trialLabel,decodeDir,times2downsample,freq,DEC,trialTypes)
    case {'eachElectrodeEachTimeEachFreq','multiElectrodeEachTimeEachFreq'}
        decode_eachTimeEachFreq_parallel(data, savefilename, trialLabel,decodeDir,times2downsample,freq,DEC,trialTypes)
    case {'eachElectrodeMultiTimeMultiFreq','multiElectrodeMultiTimeMultiFreq'}
        decode_multiTimeMultiFreq(data, savefilename, trialLabel,decodeDir,times2downsample,freq,DEC,trialTypes)
    otherwise
        keyboard
end
