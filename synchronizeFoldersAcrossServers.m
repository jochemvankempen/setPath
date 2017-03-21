% synchronize all decoding files
sub2syncLoc = {'153','154a','154b','156','168','178','180a','180b','186'};
sub2syncCFS = {'154a','154b','156','168','178'};

% sync
% 
massiveDir  = [filesep 'home' filesep 'vjochem' filesep 'Monash052' filesep 'jochem' filesep 'human' filesep 'Iowa' filesep];%localizer/data/153/TW/decode
kaniDir     = [filesep 'export' filesep 'kani' filesep 'jochem' filesep 'human' filesep 'Iowa' filesep];


% channel = sshfrommatlab('vjochem','m2.massive.org.au','1JPcoen2')


for iDataset = 2 %localizer and CFS
    
    if iDataset==1
        sub2sync = sub2syncLoc;
    elseif iDataset==2
        sub2sync = sub2syncCFS;
    end
    
    for iSub = 1:length(sub2sync)
        subID = sub2sync{iSub}
        if iDataset==1
            scp([massiveDir 'localizer' filesep 'data' filesep subID filesep 'TW' filesep 'decode' filesep '*'],[kaniDir 'localizer' filesep 'data' filesep subID filesep 'TW' filesep 'decode' filesep],1)
        elseif iDataset==2
            scptomatlab('vjochem','m2.massive.org.au','1JPcoen2',[kaniDir 'CFS' filesep 'data' filesep subID filesep 'decode' filesep],...
                [massiveDir 'CFS' filesep 'data' filesep subID filesep 'decode' filesep 'checkDecode' filesep 'TW' filesep 'high_low_visibility' filesep 'eachElectrodeEachTime' filesep 'nFoldVal_100' filesep 'ff_0_1018' filesep 'dTF_' subID '_li1_MTSPEC_high_low_visibility_cat_L1e+08.mat'])
%             syncfolder([massiveDir 'CFS' filesep 'data' filesep subID filesep 'decode' filesep ],[kaniDir 'CFS' filesep 'data' filesep subID filesep 'decode' filesep],1)
        end
    end
end

% data/153/decode/checkDecode/TW/high_low_visibility/eachElectrodeEachTime/nFoldVal_100/ff_0_1018/dTF_153_li1_MTSPEC_high_low_visibility_cat_L1e+08.mat
% channel = sshfrommatlab('vjochem','m2.massive.org.au','1JPcoen2')
channel  =  sshfrommatlabclose(channel)