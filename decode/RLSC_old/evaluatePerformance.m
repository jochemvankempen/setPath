function [P,F] =evaluatePerformance(P,F)
% 07 Mar 17 
% nClTrial( >> nClTrial(iSample,
% nTsTrial( >> nTsTrial(iSample,

pp = P;
ff = F;

PLOT.ROC = 0;

switch length(size(P))
    case 3
        ff.task = 'all';
    otherwise
end
clear P
% ff.nSession = length(ff.vSession);
% ff.nFoldValidation = length(ff.vSession);

%% evaluate performance
for iSample = 1:ff.nFoldValidation
% for iSample = 1:ff.nSession
    for iLambda = 1:length(ff.vLambda)
        switch ff.task
            case 'all'
                for iClassifier = 1:size(pp,1)
                    tmpTrain(:,iClassifier) = pp(iClassifier,iLambda,iSample).trainY';
                    tmpTest(:,iClassifier)  = pp(iClassifier,iLambda,iSample).predY';
                end

                [tmp,iTrainClass] = max(tmpTrain,[],2);
                [tmp,iTestClass] = max(tmpTest,[],2);
                if 0
                    figure(100),clf,
                    subplot(4,1,1), hold on
                    plot(tmpTrain(:,1),'ro')
                    plot(tmpTrain(:,2),'go')
                    plot(tmpTrain(:,3),'bo')
                    plot(tmpTrain(:,4),'ko')
                    axis tight
                    subplot(4,1,2), hold on
                    plot(iTrainClass,'o')
                    axis tight

                    subplot(4,1,3), hold on
                    plot(tmpTest(:,1),'ro')
                    plot(tmpTest(:,2),'go')
                    plot(tmpTest(:,3),'bo')
                    plot(tmpTest(:,4),'ko')
                    axis tight

                    subplot(4,1,4), hold on
                    plot(iTestClass,'o')
                    axis tight
                end

                for iClassifier = 1:size(pp,1)
                    vTrainTrials = ((iClassifier-1)*ff.nClTrial(iSample,1)+1) : iClassifier*ff.nClTrial(iSample,1) ;
                    corrTrain(vTrainTrials) = iTrainClass(vTrainTrials) == iClassifier ;
                    vTestTrials = ((iClassifier-1)*ff.nTsTrial(iSample,1)+1) : iClassifier*ff.nTsTrial(iSample,1) ;
                    corrTest(vTestTrials) = iTestClass(vTestTrials) == iClassifier ;
                end

                P(iLambda,iSample).iTrainClass = iTrainClass;
                P(iLambda,iSample).iTestClass = iTestClass;
                P(iLambda,iSample).corrTrain = mean(corrTrain);
                P(iLambda,iSample).corrTest = mean(corrTest);
            otherwise
                % ff.corrTrain(iLambda,iSample) = pp(iLambda,iSample).corrTrain;
                % ff.corrTest(iLambda,iSample) = pp(iLambda,iSample).corrTest;

                %% compute d'
                %%  draw an ROC curve for each lambda
                trainY = pp(iLambda,iSample).trainY;
                trainYcriterion = sort(trainY,'descend');
                clear trainPHit trainPFA
                for iCriterion = 1:length(trainYcriterion)
                    trainPHit(iCriterion) = sum(trainY(1:ff.nClTrial(iSample,1)) > trainYcriterion(iCriterion) ) /ff.nClTrial(iSample,1) ;
                    trainPFA(iCriterion)  = sum(trainY(ff.nClTrial(iSample,1)+1 : sum(ff.nClTrial(iSample,:)) ) > trainYcriterion(iCriterion) ) / ff.nClTrial(iSample,2);
                end
                ff.corrTrain(iLambda,iSample) = AreaUnderROC([trainPHit',trainPFA']);

                %% area under R
                predY = pp(iLambda,iSample).predY;
                predYcriterion = sort(predY,'descend');
                clear pHit pFA
                for iCriterion = 1:length(predYcriterion)
                    pHit(iCriterion) = sum(predY(1:ff.nTsTrial(iSample,1)) ... 
                        > predYcriterion(iCriterion) ) / ff.nTsTrial(iSample,1);
                    pFA(iCriterion) = sum(predY(ff.nTsTrial(iSample,1)+1: ... 
                        sum(ff.nTsTrial(iSample,:))) > predYcriterion(iCriterion) ) ...
                        / ff.nTsTrial(iSample,2);
                end
                ff.corrTest(iLambda,iSample) = AreaUnderROC([pHit',pFA']);

                P(iLambda,iSample).corrTest = ff.corrTest(iLambda,iSample);
                P(iLambda,iSample).corrTrain = ff.corrTrain(iLambda,iSample);
                
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

    end
end
