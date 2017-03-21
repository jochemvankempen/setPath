function [P,C] = leaveOneOutRLSC_OVA(cResp,C)
% WARNING INCOMPLETE - OBSOLETE - TO BE DOUBLE-CHECKED
%% Nov 2013
% finished, but not double checked. 
% Jochem van Kempen
%% 07 Feb 18
% one versus all classifier
%% 07 Feb 10
% Use RLSC for classification
% Compute performance for a range of lambda values.
% Decide which lambda at the end, by looking at the best test performance
% use all voxels

%% get the number of trials per session, index for each trial
% nSession = length(C.vSession);

%% compute the range of lambda

allResp = [];
for iClass = 1:length(cResp)
    allResp = [allResp;cResp{iClass}];
    nTrialPerSession(iClass) = size(cResp{iClass},1);
end

nTest = round ( min(nTrialPerSession) * C.randomSampleTest );
nTrain = min(nTrialPerSession) - nTest;

%%
% tmp = allResp*allResp';
% cLambda = round(log10(max(tmp(:))));
% vLamRange = [-5 2];
% nLam = 1;%5 ; % 15;
% % vLambda = C.vLambda;%  logspace( cLambda + vLamRange(1) , cLambda + + vLamRange(2), nLam);
% C.cLambda = cLambda;
% % C.vLamRange = vLamRange;
% C.nLam = nLam;
% % C.vLambda = vLambda;

for iSample = 1:C.nFoldValidation
    
    for iClass = 1:length(cResp)
        tmp = randperm(nTrialPerSession(iClass));
        vTestTrial{iSample,iClass} = tmp(1:nTest);
        vTrainTrial{iSample,iClass} = tmp((nTest+1):(nTest+nTrain));
        
        test {iSample,iClass} = cResp{iClass}(vTestTrial {iSample,iClass},:);
        train{iSample,iClass} = cResp{iClass}(vTrainTrial{iSample,iClass},:);
        
    end
end



for iSample = 1:C.nFoldValidation
%     if ~mod(iSample,10)
%         disp([num2str(iSample),'th sample'])
%     end
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
        
        for iLambda = 1:C.nLam
            lambda = C.vLambda(iLambda);
            w = trainRLSC(x,y,[],lambda);
            trainY = w*[ones(size(x,1),1),x]';
            predY = w*[ones(size(testX,1),1), testX]';
            
            P(iClassifier,iLambda,iSample).w = w;
            P(iClassifier,iLambda,iSample).trainY = trainY;
            P(iClassifier,iLambda,iSample).predY = predY;
        end
    end
    
end


