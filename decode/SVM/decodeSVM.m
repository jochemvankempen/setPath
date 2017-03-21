function decodeSVM(cResp,y,DP)


CC = 0.1;

% if iscell(X)
%     tmp=[];
%     for iCat = 1:length(X)
%         tmp = cat(1,tmp,X{iCat});
%     end
%     X = tmp;
% end

% make y into zeros and ones
if isempty(find(y==0))
    y(y==1)=0;
    y(y==2)=1;
end




% 07 Mar 23 Gabriel's method --
% divide data into X% (= C.randomSampleTest) for test and 100-X% for training.
% Average the performance across n repetition (Gabrield used >10)
for iClass = 1:length(cResp)
    nTrialPerSession(iClass) = size(cResp{iClass},1);
end
nTest = round ( min(nTrialPerSession) * C.randomSampleTest );
nTrain = min(nTrialPerSession) - nTest;

for iSample = 1:C.nFoldValidation
    for iClass = 1:length(cResp)
        tmp = randperm(nTrialPerSession(iClass));
        vTestTrial{iSample,iClass} = tmp(1:nTest);
        vTrainTrial{iSample,iClass} = tmp((nTest+1):(nTest+nTrain));
        
        test {iSample,iClass} = cResp{iClass}(vTestTrial {iSample,iClass},:);
        train{iSample,iClass} = cResp{iClass}(vTrainTrial{iSample,iClass},:);
        
    end
end


% create
for iSample = 1:C.nFoldValidation
    
    xTrain  = [];
    xTest   = [];
    yTrain  =[];
    yTest   =[];
    for iClass = 1:length(cResp)
        xTrain  = cat(1,xTrain,train{iSample,iClass});
        xTest   = cat(1,xTest,test{iSample,iClass});
        yTrain  = [yTrain; repmat(iClass-1,size(vTrainTrial{iSample,iClass},2),1)];
        yTest   = [yTest; repmat(iClass-1,size(vTestTrial{iSample,iClass},2),1)];
    end
    
    %train
    %     switch DP.kernelSVM
    %         case 'linear'
%     model = svmTrain(xTrain, yTrain, CC, @linearKernel);
    %         case 'gaussian'
    x1 = [1 2 1]; x2 = [0 4 -1]; sigma = 3;

    model = svmTrain(xTrain, yTrain, CC, @(x1, x2) gaussianKernel(x1, x2, sigma)); 
    %     end
    
    %test
    p{iSample} = svmPredict(model, xTest);
    pAccuracy(iSample) =  mean(double(p{iSample} == yTest)) * 100;
end

fprintf('test Accuracy: %f\n',mean(pAccuracy));
