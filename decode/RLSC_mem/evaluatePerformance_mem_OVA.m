function [P] =evaluatePerformance_mem_OVA(cResp,C,ind_train,ind_test)
% 07 Mar 17
% nClTrial( >> nClTrial(iSample,
% nTsTrial( >> nTsTrial(iSample,

% pp = P;
ff = C;

PLOT.ROC = 0;

% ff.nSession = length(ff.vSession);
% ff.nFoldValidation = length(ff.vSession);


% for iSample = 1:nSession
for iSample = 1:ff.nFoldValidation
    if mod(iSample,20)==0 & 0
        disp([num2str(iSample),'th fold validation'])
    end
    switch ff.task
        case {'FaceCheck','HappyFear','TheoryOfMind','AVsync','Varela',...
                'sfm','Gender','SacNosac','UpInv','FirstSecond','AttendNonattend','erp',...
                'fastSlowRT_CFS','CFSvsFS','cfsbm','IowaLocalizer','IowaCFS','surrogate'}
            if 0
                cl{1} = cResp{iClass}(squeeze(ind_train(iSample,iClass,1,:,1)),:);
                ts{1} = cResp{iClass}(squeeze(ind_test(iSample,iClass,1,:,1)),:);
                
                cl{2}=[];
                ts{2}=[];
                for iClassAll = classAll
                    trainIdx = squeeze(ind_train(iSample,iClass,2,:,2))==iClassAll;
                    testIdx = squeeze(ind_test(iSample,iClass,2,:,2))==iClassAll;
                    cl{2} = cat(1,cl{2},cResp{iClassAll}(squeeze(ind_train(iSample,iClass,2,trainIdx,1)),:));
                    ts{2} = cat(1,ts{2},cResp{iClassAll}(squeeze(ind_test(iSample,iClass,2,testIdx,1)),:));
                end
            elseif 0
                
                cl{1} = cResp{iClass}(squeeze(ind_train(iSample,iClass,:)),:);
                ts{1} = cResp{iClass}(squeeze(ind_test(iSample,iClass,:)),:);
                cl{2}=[];
                ts{2}=[];
                
                for iClassAll = classAll
                    
                    ts{2} = cat(1,ts{2},cResp{iClassAll}(squeeze(ind_test(iSample,iClassAll,:)),:));
                    cl{2} = cat(1,cl{2},cResp{iClassAll}(squeeze(ind_train(iSample,iClassAll,:)),:));
                end
            else
                for iClass = 1:length(cResp)
                    cl{iClass} = cResp{iClass}(squeeze(ind_train(iSample,iClass,:)),:);
                    ts{iClass} = cResp{iClass}(squeeze(ind_test(iSample,iClass,:)),:);
                end
            end
            %             cl{1} = train{iSample,1};
            %             cl{2} = train{iSample,2};
            %             ts{1} = test{iSample,1};
            %             ts{2} = test{iSample,2};
            %         case 'FaceHouse'
            %             cl{1} = [train{iSample,1};train{iSample,2}];
            %             cl{2} = [train{iSample,3};train{iSample,4}];
            %             ts{1} = [test{iSample,1};test{iSample,2}];
            %             ts{2} = [test{iSample,3};test{iSample,4}];
            %         case 'FearNeut'
            %             cl{1} = [train{iSample,1};train{iSample,3}];
            %             cl{2} = [train{iSample,2};train{iSample,4}];
            %             ts{1} = [test{iSample,1};test{iSample,3}];
            %             ts{2} = [test{iSample,2};test{iSample,4}];
            %         case 'AttendFearNeut'
            %             cl{1} = [train{iSample,1}];
            %             cl{2} = [train{iSample,2}];
            %             ts{1} = [test{iSample,1}];
            %             ts{2} = [test{iSample,2}];
            %         case 'IgnoreFearNeut'
            %             cl{1} = [train{iSample,3}];
            %             cl{2} = [train{iSample,4}];
            %             ts{1} = [test{iSample,3}];
            %             ts{2} = [test{iSample,4}];
    end
    
    [cl,ts]=normalizePatternsEachFold(cl,ts);
    
    tmpCL = cl;
    tmpTS = ts;
    clear cl ts
    for iClass = 1:length(cResp)
        classAll = find(~ismember(1:length(cResp),iClass));
        
        cl{1} = tmpCL{iClass};
        ts{1} = tmpTS{iClass};
        cl{2}=[];
        ts{2}=[];
        
        for iClassAll = classAll
            ts{2} = cat(1,ts{2},tmpTS{iClassAll});
            cl{2} = cat(1,cl{2},tmpCL{iClassAll});
        end
        
        
        [nClTrial(1)] = size(cl{1},1);
        [nClTrial(2)] = size(cl{2},1);
        [nTsTrial(1)] = size(ts{1},1);
        [nTsTrial(2)] = size(ts{2},1);
        
        ff.nClTrial = nClTrial;
        ff.nTsTrial = nTsTrial;
        
        y = [ones(nClTrial(1),1);-ones(nClTrial(2),1)];
        %    y = [ones(nClTrial(iSample,1),1);zeros(nClTrial(iSample,2),1)];
        x = [cl{1};cl{2}];
        testX = [ts{1};ts{2}];
        
        switch ff.spaceAverage
            case 1 % average of all the voxels
                x = mean(x,2);
                testX = mean(testX,2);
            case 2 % best of all the voxels
                % keyboard
                [h,p,ci,stat] = ttest2(cl{1},cl{2});
                [tmp, iVoxel] = sort(abs(stat.tstat));
                x = x(:,iVoxel(end));
                testX = testX(:,iVoxel(end));
        end
        
        
        %%
        iLambda = 1;
        % for iLambda = 1:length(ff.vLambda) % loop on lambda values has been removed
        lambda = ff.vLambda(iLambda);
        w = trainRLSC(x,y,[],lambda);
        tmpTrainY(:,iClass) = w*[ones(size(x,1),1),x]';
        tmpPredY(:,iClass) = w*[ones(size(testX,1),1), testX]';
        
        % corrTrain = ( sum(trainY(1:nClTrial(iSample,1))>0.5) +  sum(trainY(nClTrial(iSample,1)+1:sum(nClTrial))<=0.5) )/sum(nClTrial);
        % predY(1:nTsTrial(iSample,1)) = predY(1:nTsTrial(iSample,1))>0.5;
        % predY(nTsTrial(iSample,1)+1:sum(nTsTrial)) =
        % predY(nTsTrial(iSample,1)+1:sum(nTsTrial))<=0.5;
        % corrTest = mean(predY);
        if 0
            disp(['lambda =',num2str(lambda)])
            disp([num2str(corrTest) ,'/' ,num2str(corrTrain), '(test/train)'])
        end
    end
    
    [~,trainY] = max(tmpTrainY,[],2);
    [~,predY] = max(tmpPredY,[],2);
    
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
    % % for iSample = 1:ff.nSession
    %     for iLambda = 1:length(ff.vLambda)
    %         switch ff.task
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
    %% percentage correct
    for iClass = 1:iClass
        vTrainTrials = ((iClass-1)*ff.nClTrial(1)+1) : iClass*ff.nClTrial(1) ;
        tmpCorrTrain(vTrainTrials) = trainY(vTrainTrials) == iClass ;
        vTestTrials = ((iClass-1)*ff.nTsTrial(1)+1) : iClass*ff.nTsTrial(1) ;
        tmpCorrTest(vTestTrials) = predY(vTestTrials) == iClass ;
    end
    
    corrTrain(iSample) = mean(tmpCorrTrain);
    corrTest(iSample) = mean(tmpCorrTest);
    
    %             otherwise
    %                 % ff.corrTrain(iLambda,iSample) = pp(iLambda,iSample).corrTrain;
    %                 % ff.corrTest(iLambda,iSample) = pp(iLambda,iSample).corrTest;
    %
    %                 %% compute d'
    %                 %%  draw an ROC curve for each lambda
    %                 trainY = pp(iLambda,iSample).trainY;
    
    %% AUC
    %     trainYcriterion = sort(trainY,'descend');
    %     clear trainPHit trainPFA
    %     for iCriterion = 1:length(trainYcriterion)
    %         trainPHit(iCriterion) = sum(trainY(1:ff.nClTrial(1)) > trainYcriterion(iCriterion) ) /ff.nClTrial(1) ;
    %         trainPFA(iCriterion)  = sum(trainY(ff.nClTrial(1)+1 : sum(ff.nClTrial(:)) ) > trainYcriterion(iCriterion) ) / ff.nClTrial(2);
    %     end
    %     % ff.corrTrain(iLambda,iSample) = AreaUnderROC([trainPHit',trainPFA']);
    %     corrTrain(iSample) = AreaUnderROC([trainPHit',trainPFA']);
    %
    %     %% area under R
    %     % predY = pp(iLambda,iSample).predY;
    %     predYcriterion = sort(predY,'descend');
    %     clear pHit pFA
    %     for iCriterion = 1:length(predYcriterion)
    %         pHit(iCriterion) = sum(predY(1:ff.nTsTrial(1)) ...
    %             > predYcriterion(iCriterion) ) / ff.nTsTrial(1);
    %         pFA(iCriterion) = sum(predY(ff.nTsTrial(1)+1: ...
    %             sum(ff.nTsTrial(:))) > predYcriterion(iCriterion) ) ...
    %             / ff.nTsTrial(2);
    %     end
    %     % ff.corrTest(iLambda,iSample) = AreaUnderROC([pHit',pFA']);
    %     corrTest(iSample) = AreaUnderROC([pHit',pFA']);
    %
    %     %                 P(iLambda,iSample).corrTest = ff.corrTest(iLambda,iSample);
    %     %                 P(iLambda,iSample).corrTrain = ff.corrTrain(iLambda,iSample);
    %
    %     if PLOT.ROC
    %         figure(800000),clf, hold on
    %         plot(trainPFA,trainPHit,'g','linewidth',3)
    %         plot(pFA,pHit,'b','linewidth',3)
    %         plot([0 1],[0 1],'k')
    %         axis square
    %         xlabel('pFA')
    %         ylabel('pHit')
    %     end
    %%
    
end
% end
% end

%             %% evaluate
%             for iLambda = 1:ff.nLam
%                 R.corrTrain(:,iLambda) = [P(iLambda,:).corrTrain];
%                 R.corrTest(:,iLambda) = [P(iLambda,:).corrTest];
%             end

P.zmcorrTrain= mean(corrTrain);
P.zmcorrTest = mean(corrTest);
P.stdcorrTrain = std(corrTrain);
P.stdcorrTest(iClass) = std(corrTest);

%
% P.zmcorrTrain(iClass) = mean(corrTrain);
% P.zmcorrTest(iClass) = mean(corrTest);
% P.stdcorrTrain(iClass) = std(corrTrain);
% P.stdcorrTest(iClass) = std(corrTest);

