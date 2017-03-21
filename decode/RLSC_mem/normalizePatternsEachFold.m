function [cl,ts] = normalizePatternsEachFold(cl,ts)
% oldCL = cl;
% oldTS = ts;

nClass=length(cl);

% this code assumes that every class has the same number of trials
[nTrialsTr nAttributes]=size(cl{1});
[nTrialsTs nAttributes]=size(ts{1});

allPatterns=[];

for class=1:nClass
    allPatterns=[allPatterns; cl{class}];
end

meanPatterns=mean(allPatterns);
stdPatterns=std(allPatterns,0,1);

meanmat=repmat(meanPatterns,nTrialsTr,1);
stdmat=repmat(stdPatterns,nTrialsTr,1);

% nFeat = ones(size(cl{class},2),1);
% [~,zeroFeat]=find(stdmat==0);
% nFeat(zeroFeat)=0;
% nFeat = logical(nFeat);
for class=1:nClass    
%     cl{class}=(cl{class}-meanmat);
%     cl{class}(:,nFeat)=cl{class}(:,nFeat)./stdmat(:,nFeat);
    cl{class}=(cl{class}-meanmat)./stdmat;
    if sum(sum(isnan(cl{class})))>0
        error('feature NaN')
    end
end

meanmat=repmat(meanPatterns,nTrialsTs,1);
stdmat=repmat(stdPatterns,nTrialsTs,1);

% nFeat = ones(size(ts{1},2),1);
% [~,zeroFeat]=find(stdmat==0);
% nFeat(zeroFeat)=0;
% nFeat = logical(nFeat);

for class=1:nClass
%     ts{class}=(ts{class}-meanmat);
%     ts{class}(:,nFeat)=ts{class}(:,nFeat)./stdmat(:,nFeat);
    ts{class}=(ts{class}-meanmat)./stdmat;
    
    if sum(sum(isnan(ts{class})))>0
        error('feature NaN')
    end

end


% allPatternsNorm=[];
% 
% for class=1:nClass
%     allPatternsNorm=[allPatternsNorm; cl{class}];
% end
% 
% mean(allPatternsNorm)
% std(allPatternsNorm,0,1)
% 
% allPatternsTsNorm=[];
% 
% for class=1:nClass
%     allPatternsTsNorm=[allPatternsTsNorm; ts{class}];
% end
% 
% mean(allPatternsTsNorm)
% std(allPatternsTsNorm,0,1)
% 
% min(mean(allPatterns))
% max(mean(allPatterns))
% 
% min(mean(allPatternsTsNorm))
% max(mean(allPatternsTsNorm))
% 
% min(std(allPatterns,0,1))
% max(std(allPatterns,0,1))
% 
% min(std(allPatternsTsNorm,0,1))
% max(std(allPatternsTsNorm,0,1))