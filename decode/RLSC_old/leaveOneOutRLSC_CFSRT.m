function [P,C] = leaveOneOutRLSC_CFSRT(cResp,C)
% 09 Mar 25
% adopted from leaveOneOutRLSC
% use leave-one-trial-from-each-class-out procedure
% for those with many more trial number, use different strategy?

%% compute the range of lambda
% this code has been copied from the function with the same name in CFSRT/Stanford_CFS_RT_Analysis/
dbstop if error
allResp = [];
if ndims(cResp{1}) ~= 2
    disp('error : cResp should be nTrial x nFeature ')
    keyboard
end

%%
for iClass = 1:length(cResp)
    allResp = [allResp;cResp{iClass}(:)];
    [nTrialPerClass(iClass) nFeaturePerClass(iClass)  ] = size(cResp{iClass});
end
if length(unique(nFeaturePerClass)) ~= 1
    disp('error : # of features are not the same across trials')
end
if min(nTrialPerClass) < 2 ;
    disp(['# of trials are less than 2 for some class. skip.'])
    P = [];
    return
end


%%
if ~isfield(C,'nLam')
    C.nLam = 3 ;
end
if ~isfield(C,'vLamRange')
    C.vLamRange = [0 7];
    %    C.vLamRange = [-5 2];
end
if ~isfield(C,'vLambda')
    tmp = allResp*allResp';
    C.cLambda = round(log10(max(tmp(:))));
    if C.nLam == 1
        C.vLambda = logspace( C.cLambda + C.vLamRange(1) , C.cLambda + + C.vLamRange(1), C.nLam);
    else
        C.vLambda = logspace( C.cLambda + C.vLamRange(1) , C.cLambda + + C.vLamRange(2), C.nLam);
    end
else
    C.nLam = length(C.vLambda);
    % C.vLamRange = nan;
end
%%
if ~isfield(C,'nFoldValidation')
    C.nFoldValidation = min(nTrialPerClass);
end
if ~isfield(C,'randomSampleTest')
    C.randomSampleTest = []; % Gabriel used 30% for test, 70% for training.
    % this may be necessary to get better performance when the # of trials
    % are uneven (like Kanai's data, Alan's data, Hackjin's data)
end
%%
% select a leave-one-trial-for-each-class-out test and training
% create test and trial trials for nFoldValidation

for iClass = 1:length(cResp)
    % for each class, any trial can serve as a test trial only once
    vRandOrder{iClass} = randperm(nTrialPerClass(iClass));
end
%%
for iSample = 1:C.nFoldValidation
    for iClass = 1:length(cResp)
        % retrieve test trial from vRandOrder
        vTestTrial{iSample,iClass} = vRandOrder{iClass}(iSample) ;
        
        % select random training trials
        tmp = randperm(nTrialPerClass(iClass));
        tmp = tmp(1:min(nTrialPerClass)); % select first 1:nTrial as training
        if sum(ismember(tmp,vTestTrial{iSample,iClass})) == 0
            % training did not contain the test trial
            tmp = tmp(1:end-1);
        else
            % training contained the test trial
            tmp = tmp(~ismember(tmp,vTestTrial{iSample,iClass}));
        end
        if length(union(tmp,vTestTrial{iSample,iClass})) ~= min(nTrialPerClass)
            disp('error')
            keyboard
        end
        vTrainTrial{iSample,iClass} = tmp;
        
        test {iSample,iClass} = cResp{iClass}(vTestTrial{iSample,iClass},:);
        train{iSample,iClass} = cResp{iClass}(vTrainTrial{iSample,iClass},:);
        
    end
end

%%
nTest = 1;
nTrain = size(train{1,1},1);

% keyboard;

%%
for iSample = 1:C.nFoldValidation
    disp([num2str(iSample),'th fold validation'])    
%     if mod(iSample,20)==0 & 0
%         disp([num2str(iSample),'th fold validation'])
%     end
    
    x = [];
    testX = [];
    for iClass = 1:length(cResp)
        x = [x; train{iSample,iClass}];
        testX = [testX; test{iSample,iClass}];
        
        [nClTrial(iSample,iClass)] = size(train{iSample,iClass},1);
        [nTsTrial(iSample,iClass)] = size(test{iSample,iClass},1);
    end
    
    C.nClTrial = nClTrial;
    C.nTsTrial = nTsTrial;
    
    for iClassifier = 1:length(cResp)
        y = [];
        
        for iCl = 1:length(cResp)
            if iCl == iClassifier
                y = [y;ones(nTrain,1);];
            else
                y = [y;zeros(nTrain,1);];
            end
        end
        
        for iLambda = 1:length(C.vLambda)
            lambda = C.vLambda(iLambda);
            w = trainRLSC(x,y,[],lambda);
            trainY = w*[ones(size(x,1),1),x]';
            predY = w*[ones(size(testX,1),1), testX]';
            
            P(iLambda,iSample).w(iClassifier,:) = w;
            P(iLambda,iSample).trainY(iClassifier,:) = trainY;
            P(iLambda,iSample).predY(iClassifier,:) = predY;
        end
    end
    
end


