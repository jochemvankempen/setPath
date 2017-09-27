function [P] =evaluatePerformance_mem(cResp,DEC,ind_train,ind_test)
% 07 Mar 17
% nClTrial( >> nClTrial(iSample,
% nTsTrial( >> nTsTrial(iSample,

% pp = P;

PLOT.ROC = 0;

% ff.nSession = length(ff.vSession);
% ff.nFoldValidation = length(ff.vSession);


% for iSample = 1:nSession
for iSample = 1:DEC.nCrossVal
    if mod(iSample,20)==0 & 0
        disp([num2str(iSample),'th fold validation'])
    end
    
    cl = [];
    cl{1} = cResp{1}(squeeze(ind_train(iSample,1,:)),:);
    cl{2} = cResp{2}(squeeze(ind_train(iSample,2,:)),:);
    ts{1} = cResp{1}(squeeze(ind_test(iSample,1,:)),:);
    ts{2} = cResp{2}(squeeze(ind_test(iSample,2,:)),:);
    
    [cl,ts]=normalizePatternsEachFold(cl,ts);
    
    [nClTrial(1)] = size(cl{1},1);
    [nClTrial(2)] = size(cl{2},1);
    [nTsTrial(1)] = size(ts{1},1);
    [nTsTrial(2)] = size(ts{2},1);
    
    DEC.nClTrial = nClTrial;
    DEC.nTsTrial = nTsTrial;
    
    y = [ones(nClTrial(1),1);-ones(nClTrial(2),1)];
    %    y = [ones(nClTrial(iSample,1),1);zeros(nClTrial(iSample,2),1)];
    x = [cl{1};cl{2}];
    testX = [ts{1};ts{2}];
        
    
    %%
    w = trainRLSC(x,y,[],DEC.lambda);
    trainY = w*[ones(size(x,1),1),x]';
    predY = w*[ones(size(testX,1),1), testX]';
    
    % corrTrain = ( sum(trainY(1:nClTrial(iSample,1))>0.5) +  sum(trainY(nClTrial(iSample,1)+1:sum(nClTrial))<=0.5) )/sum(nClTrial);
    % predY(1:nTsTrial(iSample,1)) = predY(1:nTsTrial(iSample,1))>0.5;
    % predY(nTsTrial(iSample,1)+1:sum(nTsTrial)) =
    % predY(nTsTrial(iSample,1)+1:sum(nTsTrial))<=0.5;
    % corrTest = mean(predY);
    if 0
        disp(['lambda =',num2str(lambda)])
        disp([num2str(corrTest) ,'/' ,num2str(corrTrain), '(test/train)'])
    end
    
    %         P(iLambda,iSample).w = w;
    %         P(iLambda,iSample).trainY = trainY;
    %         P(iLambda,iSample).predY = predY;
    
    
    
    % P(iLambda,iSample).corrTrain = corrTrain;
    % P(iLambda,iSample).corrTest = corrTest;
    
    %         if corrTrain <= corrTest & corrTrain < 1
    %             disp('too much penalty')
    %             break
    %         end
    
    
    % end
    
    
    
    
    
    
    
    
    
    % %% evaluate performance
    % % for iSample = 1:DEC.nSession
    %     for iLambda = 1:length(DEC.vLambda)
    %         switch DEC.task
    %             case 'all'
    %                 for iClassifier = 1:4
    %                     tmpTrain(:,iClassifier) = pp(iClassifier,iLambda,iSample).trainY';
    %                     tmpTest(:,iClassifier)  = pp(iClassifier,iLambda,iSample).predY';
    %                 end
    %
    %                 [tmp,iTrainClass] = max(tmpTrain,[],2);
    %                 [tmp,iTestClass] = max(tmpTest,[],2);
    %                 if 0
    %                     figure(100),clf,
    %                     subplot(4,1,1), hold on
    %                     plot(tmpTrain(:,1),'ro')
    %                     plot(tmpTrain(:,2),'go')
    %                     plot(tmpTrain(:,3),'bo')
    %                     plot(tmpTrain(:,4),'ko')
    %                     axis tight
    %                     subplot(4,1,2), hold on
    %                     plot(iTrainClass,'o')
    %                     axis tight
    %
    %                     subplot(4,1,3), hold on
    %                     plot(tmpTest(:,1),'ro')
    %                     plot(tmpTest(:,2),'go')
    %                     plot(tmpTest(:,3),'bo')
    %                     plot(tmpTest(:,4),'ko')
    %                     axis tight
    %
    %                     subplot(4,1,4), hold on
    %                     plot(iTestClass,'o')
    %                     axis tight
    %                 end
    %
    %                 for iClassifier = 1:4
    %                     vTrainTrials = ((iClassifier-1)*DEC.nClTrial(iSample,1)+1) : iClassifier*DEC.nClTrial(iSample,1) ;
    %                     corrTrain(vTrainTrials) = iTrainClass(vTrainTrials) == iClassifier ;
    %                     vTestTrials = ((iClassifier-1)*DEC.nTsTrial(iSample,1)+1) : iClassifier*DEC.nTsTrial(iSample,1) ;
    %                     corrTest(vTestTrials) = iTestClass(vTestTrials) == iClassifier ;
    %                 end
    %
    %                 P(iLambda,iSample).iTrainClass = iTrainClass;
    %                 P(iLambda,iSample).iTestClass = iTestClass;
    %                 P(iLambda,iSample).corrTrain = mean(corrTrain);
    %                 P(iLambda,iSample).corrTest = mean(corrTest);
    %             otherwise
    %                 % DEC.corrTrain(iLambda,iSample) = pp(iLambda,iSample).corrTrain;
    %                 % DEC.corrTest(iLambda,iSample) = pp(iLambda,iSample).corrTest;
    %
    %                 %% compute d'
    %                 %%  draw an ROC curve for each lambda
    %                 trainY = pp(iLambda,iSample).trainY;
    trainYcriterion = sort(trainY,'descend');
    clear trainPHit trainPFA
    for iCriterion = 1:length(trainYcriterion)
        trainPHit(iCriterion) = sum(trainY(1:DEC.nClTrial(1)) > trainYcriterion(iCriterion) ) /DEC.nClTrial(1) ;
        trainPFA(iCriterion)  = sum(trainY(DEC.nClTrial(1)+1 : sum(DEC.nClTrial(:)) ) > trainYcriterion(iCriterion) ) / DEC.nClTrial(2);
    end
    % DEC.corrTrain(iLambda,iSample) = AreaUnderROC([trainPHit',trainPFA']);
    corrTrain(iSample) = AreaUnderROC([trainPHit',trainPFA']);
    %% area under R
    % predY = pp(iLambda,iSample).predY;
    predYcriterion = sort(predY,'descend');
    clear pHit pFA
    for iCriterion = 1:length(predYcriterion)
        pHit(iCriterion) = sum(predY(1:DEC.nTsTrial(1)) ...
            > predYcriterion(iCriterion) ) / DEC.nTsTrial(1);
        pFA(iCriterion) = sum(predY(DEC.nTsTrial(1)+1: ...
            sum(DEC.nTsTrial(:))) > predYcriterion(iCriterion) ) ...
            / DEC.nTsTrial(2);
    end
    % DEC.corrTest(iLambda,iSample) = AreaUnderROC([pHit',pFA']);
    corrTest(iSample) = AreaUnderROC([pHit',pFA']);
    
    %                 P(iLambda,iSample).corrTest = DEC.corrTest(iLambda,iSample);
    %                 P(iLambda,iSample).corrTrain = DEC.corrTrain(iLambda,iSample);
    
    if PLOT.ROC
        figure(800000),clf, hold on
        plot(trainPFA,trainPHit,'g','linewidth',3)
        plot(pFA,pHit,'b','linewidth',3)
        plot([0 1],[0 1],'k')
        axis square
        xlabel('pFA')
        ylabel('pHit')
    end
    
end

% end
% end

%             %% evaluate
%             for iLambda = 1:DEC.nLam
%                 R.corrTrain(:,iLambda) = [P(iLambda,:).corrTrain];
%                 R.corrTest(:,iLambda) = [P(iLambda,:).corrTest];
%             end

P.zmcorrTrain = mean(corrTrain);
P.zmcorrTest = mean(corrTest);
P.stdcorrTrain = std(corrTrain);
P.stdcorrTest = std(corrTest);