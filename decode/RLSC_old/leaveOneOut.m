function [TS, TR] = leaveOneOut(allResp,C)
[nTrial(1) nVoxel] = size(allResp{1});
[nTrial(2) nVoxel] = size(allResp{2});
nAllTrial = sum(nTrial);
nRmAllTrial = nAllTrial -1;
nUsedVoxel = C.nUsedVoxel;

nSession = length(C.vSession);

for iVoxel = C.vUsedVoxel
    % compute covariance matrix, by leaving out one sample
    disp(['# of voxels used = ',num2str(iVoxel)])
    for iSample = 1:sum(nTrial(1:2))
        if mod(iSample,20)==0
            disp([num2str(iSample),'th sample'])
        end

        % select a leave-one-out sample and remaining two classes      
        if iSample <= nTrial(1)
            sample = allResp{1}(iSample,:);
            cl{1} = [allResp{1}(setdiff([1:nTrial(1)],iSample),:);] ;
            cl{2} = [allResp{2}];
        else
            iTmpSample = iSample - nTrial(1);
            sample = allResp{2}(iTmpSample,:);
            cl{1} = [allResp{1}];
            cl{2} = [allResp{2}(setdiff([1:nTrial(2)],iTmpSample),:)];
        end

        [nClTrial(1)] = size(cl{1},1);
        [nClTrial(2)] = size(cl{2},1);

        % get the tscore for discrimination -- >> use the most
        % discriminative iVoxels
        [h,p,ci,stats] = ttest2(cl{1},cl{2});
        [sorted_tstat, i_sort_tstat] = sort(abs(stats.tstat),2,'descend');
        % [sorted_tstat, i_sort_tstat] = sort(stats.tstat,2,'ascend');
        % [sorted_tstat, i_sort_tstat] = sort(stats.tstat,2,'descend');
        vUseVoxels = i_sort_tstat(1:iVoxel);

        %% 06 Dec 26 -- fixed. Before pooling, subtract the means!
        % create pooled training data
        % rmAll = [cl{1}(:,vUseVoxels);cl{2}(:,vUseVoxels)];

        for iClass = 1:2
            tmpcl{iClass} = cl{iClass}(:,vUseVoxels)-repmat(mean(cl{iClass}(:,vUseVoxels)),[nClTrial(iClass) 1]);
        end
        rmAll = [tmpcl{1};tmpcl{2}];
        pooledCov = cov(rmAll);
        invCov = inv(pooledCov);

        % calculate exponent y for the posterior probability >> see
        % Hampton & ODoherty (in Press)
        for iClass = 1:2
            TS.y(iVoxel,iSample,iClass) = -1/2*(sample(:,vUseVoxels) - mean(cl{iClass}(:,vUseVoxels)))*invCov*(sample(:,vUseVoxels) - mean(cl{iClass}(:,vUseVoxels)))'+log(nClTrial(1)/nRmAllTrial);
        end

        % decide which class to assign.  use discriminant function
        % -- can be made faster by calculate this as matrix operation
        if iSample <= nTrial(1)
            TS.yy(iVoxel,iSample)    = TS.y(iVoxel,iSample,1)-TS.y(iVoxel,iSample,2);
        else
            TS.yy(iVoxel,iSample)    = TS.y(iVoxel,iSample,2)-TS.y(iVoxel,iSample,1);
        end


        if C.calculateTrainingError
            iTrSampleAll = 0;
            for iTrClass = 1:length(cl)
                for iTrSample = 1:nClTrial(iTrClass)
                    iTrSampleAll = iTrSampleAll + 1;
                    trSample = cl{iTrClass}(iTrSample,vUseVoxels);

                    % TS.y for training sample
                    for iClass = 1:2
                        TR.y(iVoxel,iTrSampleAll,iSample,iClass) = -1/2*(trSample - mean(cl{iClass}(:,vUseVoxels)))*invCov*(trSample - mean(cl{iClass}(:,vUseVoxels)))'+log(nClTrial(iTrClass)/nRmAllTrial);
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
                    [h,p,ci,allstats] = ttest2(allResp{1},allResp{2});
                    [allsorted_tstat, i_sort_tstat] = sort((allstats.tstat),2,'descend');
                    plot(allsorted_tstat)
                    title(['tscore: mean=',num2str(mean(allstats.tstat),2) ': median=',num2str(median(allstats.tstat),2)])
                end
            end

            if iVoxel == 2 & iSample == 1
                subplot(2,2,2)
                hold on
                plot(cl{1}(:,vUseVoxels(1)),cl{1}(:,vUseVoxels(2)),'o')
                plot(cl{2}(:,vUseVoxels(1)),cl{2}(:,vUseVoxels(2)),'rx')
                plot(sample(:,vUseVoxels(1)),sample(:,vUseVoxels(2)),'g*');

                X = [sample(:,vUseVoxels);cl{1}(:,vUseVoxels);cl{2}(:,vUseVoxels)]';
                [D,N] = size(X);

                % STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS
                X2 = sum(X.^2,1);
                distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
                [sorted,index] = sort(distance);

                subplot(2,2,3),cla
                n1 = nClTrial(1)+1;
                imagesc(distance), hold on
                plot([n1,n1]+.5,[0 N],'w')
                plot([0 N],[n1,n1]+.5,'w')
            end
        end

    end
    TS.pCorr1(iVoxel) = mean(TS.yy(iVoxel,1:nTrial(1))>0);
    TS.pCorr2(iVoxel) = mean(TS.yy(iVoxel,nTrial(1)+1:end)>0);
    TS.pCorr(iVoxel)  = mean([TS.pCorr1(iVoxel),TS.pCorr2(iVoxel)]);

    disp(['correct test/training classification = ',num2str(TS.pCorr(iVoxel)) '/' num2str(mean(TR.pCorr(iVoxel,:)))])
end

% used for averageClassificationPerformance.m
TS.nTrial = nTrial;
