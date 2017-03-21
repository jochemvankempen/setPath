function cResp = normalizePatterns(cResp)

nClass=length(cResp);


allPatterns=[];

for class=1:nClass
    allPatterns=[allPatterns; cResp{class}];
end

meanPatterns=mean(allPatterns);
stdPatterns=std(allPatterns,0,1);

for class=1:nClass
    [nTrials nAttributes]=size(cResp{class});
    
    meanmat=repmat(meanPatterns,nTrials,1);
    stdmat=repmat(stdPatterns,nTrials,1);
    
    cResp{class}=(cResp{class}-meanmat)./stdmat;
end