function [TS, TR] = leaveOneOutSession(cResp,C)
nSession = length(C.vSession);
%% get the number of trials per session, index for each trial
% % vTrials=[];% for index
% vSubTrials = [];
% for iSession = 1:nSession
%     tmp_vTrials = [1:12]+(iSession-1)*12*4;
%     tmp_vTrials = repmat(tmp_vTrials,[4 1]) + repmat([0 12 24 36]',[1,12]);
%     vSubTrials = [vSubTrials,repmat([1:12],[4 1])+repmat([0 12 24 36]',[1,12])];
% %    vTrials = [vTrials,tmp_vTrials];
% end
%%

for iVoxel = C.vUsedVoxel
    % compute covariance matrix, by leaving out one sample
    disp(['# of voxels used = ',num2str(iVoxel)])
    %    for iSample = 1:sum(nTrial(1:2))
    for iSample = 1:nSession
        if mod(iSample,20)==0
            disp([num2str(iSample),'th sample'])
        end

        % select a leave-one-out test and training
        vTrial = [1:12] +(iSample-1)*12;
        for iClass = 1:4
            test{iClass}= cResp{iClass}(vTrial,:);
            train{iClass} = cResp{iClass}(setdiff([1:size(cResp{iClass},1)],vTrial),:);
        end

        cl = [];
        switch C.task
            case 'FaceHouse'
                cl{1} = [train{1};train{2}];
                cl{2} = [train{3};train{4}];
                ts{1} = [test{1};test{2}];
                ts{2} = [test{3};test{4}];
            case 'FearNeut'
                cl{1} = [train{1};train{3}];
                cl{2} = [train{2};train{4}];
                ts{1} = [test{1};test{3}];
                ts{2} = [test{2};test{4}];
            case 'AttendFearNeut'
                cl{1} = [train{1}];
                cl{2} = [train{2}];
            case 'IgnoreFearNeut'
                cl{1} = [train{3}];
                cl{2} = [train{4}];
        end

        [nClTrial(1)] = size(cl{1},1);
        [nClTrial(2)] = size(cl{2},1);
        [nTsTrial(1)] = size(ts{1},1);
        [nTsTrial(2)] = size(ts{2},1);





        switch C.classifier
            case 'RLSC'

                
                y = [ones(nClTrial(1),1);zeros(nClTrial(2),1)];
                x = [cl{1};cl{2}];
                                    testX = [ts{1};ts{2}];

                
                tmp = x*x';
                cLambda = round(log10(max(tmp(:))));
                vLamRange = [-5 2]
                nLam = 15
                vLambda = logspace( cLambda + vLamRange(1) , cLambda + + vLamRange(2), nLam); 
                
                for iLambda = 1:length(vLambda)
                    lambda = vLambda(iLambda);
                    w(iLambda,:) = trainRLSC(x,y,[],lambda);
                    trainY = w(iLambda,:)*[ones(size(x,1),1),x]';
                    corrTrain(iLambda) = ( sum(trainY(1:nClTrial(1))>0.5) +  sum(trainY(nClTrial(1)+1:sum(nClTrial))<=0.5) )/sum(nClTrial);

                    predY = w(iLambda,:)*[ones(size(testX,1),1), testX]';
                    predY(1:nTsTrial(1)) = predY(1:nTsTrial(1))>0.5;
                    predY(nTsTrial(1)+1:sum(nTsTrial)) = predY(nTsTrial(1)+1:sum(nTsTrial))<=0.5;

                    corrTest(iLambda) = mean(predY);
                    disp(['lambda =',num2str(lambda)])
                    %        disp(['weight =' num2str(w)])
                    disp([num2str(corrTest(iLambda)) ,'/' ,num2str(corrTrain(iLambda)), '(test/train)'])

                end

            case 'mahal'
                % get the tscore for discrimination -- >> use the most
                % discriminative iVoxels
                [h,p,ci,stats] = ttest2(cl{1},cl{2});
                [sorted_tstat, i_sort_tstat] = sort(abs(stats.tstat),2,'descend');
                % [sorted_tstat, i_sort_tstat] = sort(stats.tstat,2,'ascend');
                % [sorted_tstat, i_sort_tstat] = sort(stats.tstat,2,'descend');
                vUseVoxels = i_sort_tstat(1:iVoxel);
                switch C.task
                    case 'FaceHouse'
                        ts{1} = [test{1}(:,vUseVoxels);test{2}(:,vUseVoxels)];
                        ts{2} = [test{3}(:,vUseVoxels);test{4}(:,vUseVoxels)];
                    case 'FearNeut'
                        ts{1} = [test{1}(:,vUseVoxels);test{3}(:,vUseVoxels)];
                        ts{2} = [test{2}(:,vUseVoxels);test{4}(:,vUseVoxels)];
                    case 'AttendFearNeut'
                        ts{1} = [test{1}(:,vUseVoxels)];
                        ts{2} = [test{2}(:,vUseVoxels)];
                    case 'IgnoreFearNeut'
                        ts{1} = [test{3}(:,vUseVoxels)];
                        ts{2} = [test{4}(:,vUseVoxels)];
                end
                if 0
                    figure
                    subplot(2,1,1)
                    hist([cl{1}(:,vUseVoxels);cl{2}(:,vUseVoxels)])
                    subplot(2,1,2)
                    hist([ts{1};ts{2}])
                end
                for iClass = 1:length(cl)
                    tr{iClass} = cl{iClass}(:,vUseVoxels);
                end
                %% 06 Dec 26 -- fixed. Before pooling, subtract the means!
                % create pooled training data
                % rmAll = [cl{1}(:,vUseVoxels);cl{2}(:,vUseVoxels)];
                for iClass = 1:2
                    tmpcl{iClass} = tr{iClass} -repmat(mean(tr{iClass}),[nClTrial(iClass) 1]);
                end
                rmAll = [tmpcl{1};tmpcl{2}];
                pooledCov = cov(rmAll);
                invCov = inv(pooledCov);


                % decide which class to assign.  use discriminant function
                % -- can be made faster by calculate this as matrix operation
                %         if iSample <= nTrial(1)
                %             TS.yy(iVoxel,iSample)    = TS.y(iVoxel,iSample,1)-TS.y(iVoxel,iSample,2);
                %         else
                %             TS.yy(iVoxel,iSample)    = TS.y(iVoxel,iSample,2)-TS.y(iVoxel,iSample,1);
                %         end

                % TS.yy(iVoxel,iSample,1
                iTestClass = 1;
                tmpyy(iTestClass,:) = (tmpy(1,iTestClass,:)-tmpy(2,iTestClass,:)); % for iTestClass == 1, distance to mean of iTrainingClass1 should be smaller
                iTestClass = 2;
                tmpyy(iTestClass,:) = (tmpy(2,iTestClass,:)-tmpy(1,iTestClass,:)); % for iTestClass == 2, distance to mean of iTrainingClass2 should be smaller

                pCorr(iVoxel,iSample,:) = mean(tmpyy'>0) ;

                TS.yy(iVoxel,iSample,:,:) = tmpyy; % the third demension is for test class, forth demension is for # of test trials

                if C.calculateTrainingError
                    iTrSampleAll = 0;
                    for iTrClass = 1:length(cl)
                        for iTrSample = 1:nClTrial(iTrClass)
                            iTrSampleAll = iTrSampleAll + 1;
                            trSample = cl{iTrClass}(iTrSample,vUseVoxels);

                            % TS.y for training sample
                            for iClass = 1:2
                                tmp = trSample - mean(tr{iClass});
                                TR.y(iVoxel,iTrSampleAll,iSample,iClass) = -1/2*tmp*invCov*tmp'+log(nClTrial(iTrClass)/sum(nClTrial));
                            end
                            if iTrClass == 1
                                TR.yy(iVoxel,iTrSampleAll,iSample) = TR.y(iVoxel,iTrSampleAll,iSample,1)-TR.y(iVoxel,iTrSampleAll,iSample,2);
                            else
                                TR.yy(iVoxel,iTrSampleAll,iSample) = TR.y(iVoxel,iTrSampleAll,iSample,2)-TR.y(iVoxel,iTrSampleAll,iSample,1);
                            end
                        end
                    end

                    % % of correct is the average for correct classification in
                    % class 1 and class 2.  >> see Hampton & ODoherty
                    TR.pCorr1(iVoxel,iSample) = mean(TR.yy(iVoxel,1:nClTrial(1),iSample)>0);
                    TR.pCorr2(iVoxel,iSample) = mean(TR.yy(iVoxel,nClTrial(1)+1:end,iSample)>0);
                    TR.pCorr(iVoxel,iSample) = mean([TR.pCorr1(iVoxel,iSample),TR.pCorr2(iVoxel,iSample)]);

                    TR.ppCorr1(iVoxel,iSample) = 1./(1+exp(-mean(TR.yy(iVoxel,1:nClTrial(1),iSample))));
                    TR.ppCorr2(iVoxel,iSample) = 1./(1+exp(-mean(TR.yy(iVoxel,nClTrial(1)+1:end,iSample))));
                    TR.ppCorr(iVoxel,iSample) = mean([TR.ppCorr1(iVoxel,iSample),TR.ppCorr2(iVoxel,iSample)]);

                end

                if C.plotTscoreDistance
                    if  iVoxel == 1 & iSample == 1
                        subplot(2,2,1)
                        hold on
                        if 0
                            plot(cl{1}(:,vUseVoxels),'o')
                            plot(cl{2}(:,vUseVoxels),'rx')
                            plot(sample(:,vUseVoxels),'g*');
                        else
                            [h,p,ci,allstats] = ttest2([cResp{1},cResp{2}]',[cResp{3},cResp{4}]');
                            [allsorted_tstat, i_sort_tstat] = sort((allstats.tstat),2,'descend');
                            plot(allsorted_tstat)
                            title(['tscore: mean=',num2str(mean(allstats.tstat),2) ': median=',num2str(median(allstats.tstat),2)])
                        end
                    end

                    if iVoxel == 2 & iSample == 1
                        subplot(2,2,2)
                        hold on
                        plot(cl{1}(:,vUseVoxels(1)),cl{1}(:,vUseVoxels(2)),'bo')
                        plot(cl{2}(:,vUseVoxels(1)),cl{2}(:,vUseVoxels(2)),'rx')
                        plot(ts{1}(:,1),ts{1}(:,2),'g*');
                        plot(ts{2}(:,1),ts{2}(:,2),'ys');

                        X = [ts{1};ts{2};cl{1}(:,vUseVoxels);cl{2}(:,vUseVoxels)]';
                        [D,N] = size(X);

                        % STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS
                        X2 = sum(X.^2,1);
                        distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
                        [sorted,index] = sort(distance);

                        subplot(2,2,3),cla
                        n1 = [nTest, nTest*2 nTest*2+nClTrial(1)+1];
                        imagesc(distance), hold on
                        plot([n1;n1]+.5,repmat([0 N],[3 1])','w')
                        plot(repmat([0 N],[3 1])',[n1;n1]+.5,'w')
                    end
                end
            TS.pCorr(iVoxel,:,:)  = squeeze(mean(TS.yy(iVoxel,:,:,:)>0,4));

        end
    end

    %     TS.pCorr1(iVoxel) = mean(TS.yy(iVoxel,1:nTrial(1))>0);
    %     TS.pCorr2(iVoxel) = mean(TS.yy(iVoxel,nTrial(1)+1:end)>0);
    %     TS.pCorr(iVoxel)  = mean([TS.pCorr1(iVoxel),TS.pCorr2(iVoxel)]);

    disp(['correct test/training classification = ',num2str(mean(mean(TS.pCorr(iVoxel,:,:)))), '/' num2str(mean(TR.pCorr(iVoxel,:)))])
end

% used for averageClassificationPerformance.m
% TS.nTrial = nTrial;
