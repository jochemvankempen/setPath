%% leave one out validation
% svm
% strange why the performance doesn't drop to .5 for 1 voxel
% n-fold validation
% create test and training at once

clear all
dbstop if error
tic
rand('state',0)

RLSC.compute = 1;
RLSC.sortByTtest= 0;
RLSC.calculateTrainingError= 0;
RLSC.printFigure = 1;

SVM = 0;
PLOT.RLSC = 1;

PLOT.covariance = 1;
PLOT.ttest = 0;

nVoxel = 1000;
separation = .1;
mu{1} = separation*randn(nVoxel,1).*ones(nVoxel,1);
mu{2} = -separation*randn(nVoxel,1).*ones(nVoxel,1);
correlation = 0.1;
% S = correlation*ones(nVoxel)+eye(nVoxel)*(1-correlation);
S = correlation*rand(nVoxel)+eye(nVoxel);
S = S'*S;
S = S/max(S(:));

if PLOT.covariance
    figure(100),
    imagesc(S),colorbar
end

nData = 320;
nFoldValidation = 8;
nTest = nData/nFoldValidation
nTrain = nData - nTest

colors = 'rgbk'
markers = 'oxox'
OUTLIER = 0 % 5;

omu{1} = zeros(nVoxel,1);
omu{1}(1) = 4;
omu{2} = zeros(nVoxel,1);
omu{2}(1) = -4;
os = S;

%%
figure(1),clf, hold on
for i = 1:2
    data{i} = mvnrnd(mu{i},S,nData);
    if OUTLIER
        data{i}(1:OUTLIER,:) = mvnrnd(omu{i},S,OUTLIER);
    end

    plot(data{i}(:,1),data{i}(:,2),'color',colors(i),'marker',markers(i),'linestyle','none')
end
colors = 'bk';
markers = 'ox';
axis tight
vx = xlim
vy = ylim

%%

PLOT.boundary = 0;
saveW = 0 ;
clear corrTest corrTrain w
scale = 'log'
nLambda = 1
switch scale
    case 'linear'
        vLambda = linspace(0,1,nLambda)
    case 'log'
        vLambda = [0 logspace(-nLambda+1,0,nLambda)]
end

% keyboard
% vUsedVoxel = 1:10:100% nVoxel
%%

% vUsedVoxel = 1:10:nVoxel
if RLSC.sortByTtest
%     vUsedVoxel = [1:40] %
    vUsedVoxel = [1:10:100 200 400 1000] %
else
    vUsedVoxel = [1:10:100 200 400 1000] %
    ori_varsToUse = randperm(nVoxel);
end
% 1:200:nVoxel

%% validation

vTrial = 1:nData;


if RLSC.compute
    for iValidation = 1:nFoldValidation
        vTest = [1:nData/nFoldValidation] + (iValidation-1)*nData/nFoldValidation;
        vTrain =vTrial(~ismember(vTrial,vTest));
        y = [ones(nTrain,1);zeros(nTrain,1)];

        if length(unique([vTest, vTrain])) ~= nData
            error('')
        end

        for iClass = 1:2
            train{iClass} = data{iClass}(vTrain,:);
            test{iClass} = data{iClass}(vTest,:);
        end

        for iLambda = 1:length(vLambda)
            lambda = vLambda(iLambda)
            datestr(now)

            if RLSC.sortByTtest
                [h,p,ci,stat] = ttest(train{1},train{2});
                [sortTstat,iSortTstat] = sort(abs(stat.tstat));
            end
            for iVoxel = vUsedVoxel
                if RLSC.sortByTtest
                    varsToUse = iSortTstat(1:iVoxel);
                else
                    varsToUse = ori_varsToUse(1:iVoxel);
                end
                x = [train{1};train{2}];


                w = trainRLSC(x,y,varsToUse,lambda);
                if PLOT.boundary
                    % clear h
                    h(1) = plot(vx,-vx*w(2)/w(3)-w(1)/w(3),'g-','linewidth',iLambda);
                    h(2) = plot(vx,-vx*w(2)/w(3)-w(1)/w(3)+1/w(3),'r-','linewidth',iLambda);
                    h(3) = plot(vx,-vx*w(2)/w(3)-w(1)/w(3)+.5/w(3),'b-','linewidth',iLambda);
                    % pause
                end
                %%
                %%
                % test error
                if RLSC.calculateTrainingError
                    trainY = w*[ones(size(x,1),1), x(:,varsToUse)]';
                    corrTrain(iLambda,iVoxel,iValidation) = ( sum(trainY(1:nTrain)>0.5) +  sum(trainY(nTrain+1:2*nTrain)<=0.5) )/nTrain/2;
                end
                %    disp(corrTrain/nTrain/2)
                %%
                testX = [test{1};test{2}];
                testY = w*[ones(size(testX,1),1), testX(:,varsToUse)]';
                corrTest(iLambda,iVoxel,iValidation) = ( sum(testY(1:nTest)>0.5) +  sum(testY(nTest+1:2*nTest)<=0.5) ) /nTest/2 ;
                %    disp(corrTest/nTest/2)
                if RLSC.calculateTrainingError
                    disp([num2str(iVoxel),' voxels: ' ,num2str(corrTest(iLambda,iVoxel,iValidation)) ,'/' ,num2str(corrTrain(iLambda,iVoxel,iValidation)), '(test/train)'])
                else
                    disp([num2str(iVoxel),' voxels: ' ,num2str(corrTest(iLambda,iVoxel,iValidation)) ,'(test)'])
                end


                if saveW
                    W{iVoxel}(iLambda,:) = w;
                end
            end
        end
    end
end
%%
if PLOT.RLSC
    figure,clf,hold on,
    colors = 'rgbkymc'

    meanCorrTest = mean(corrTest,3);
    stderrCorrTest = std(corrTest,0,3) / sqrt(size(corrTest,3));
    
    for iLambda = 1:length(vLambda)
        plot(vUsedVoxel,meanCorrTest(iLambda,vUsedVoxel),'o-','color',colors(iLambda))
        legendtext{iLambda} = ['test lm=',num2str(vLambda(iLambda))];
    end

    ms = meanCorrTest - stderrCorrTest;
    ps = meanCorrTest + stderrCorrTest;
    
    for iLambda = 1:length(vLambda)
        plot([vUsedVoxel;vUsedVoxel],[ms(iLambda,vUsedVoxel); ps(iLambda,vUsedVoxel)],'-','color',colors(iLambda))
    end
            
    
    if RLSC.calculateTrainingError
        for iLambda = 1:length(vLambda)
            plot(vUsedVoxel,mean(corrTrain(iLambda,vUsedVoxel,:),3),'x:','color',colors(iLambda))
            legendtext{iLambda+length(vLambda)} = ['train lm=',num2str(vLambda(iLambda))];
        end
    end
    legend(legendtext,'location','northwest')
    
   ylabel('% corr. classification')
   xlabel('# of voxles used')
    
    %plot(vLambda,corrTrain,'bo-')
    %plot(vLambda,corrTest,'kx-')
    if RLSC.sortByTtest
%        set(gca,'xscale','log')
        title('voxels selected by ttest')
    else
        title('voxels chosen at random')        
    end
    
    if RLSC.printFigure
        print(gcf,'-djpeg',['simulationFMRI',num2str(RLSC.sortByTtest),'.jpg'])
    end
end
%%
if PLOT.ttest
    [h,p,ci,stat] = ttest(data{1},data{2});
    %%
    [h,p,ci,teststat] = ttest(test{1},test{2});
    %%
    figure(4),clf
    [ sortedTrainTstat,iSortedTrainTstat ] = sort(stat.tstat);
    plot(sortedTrainTstat)
    %%
    figure(5)
    plot(teststat.tstat(iSortedTrainTstat))
end


%% svm
if SVM
    varsToUse = ori_varsToUse
    gamma = logspace(-15,5,10)
    cost = logspace(-5,15,10)

    y = [ones(nData,1);zeros(nData,1)];
    x = [data{1}(:,varsToUse);data{2}(:,varsToUse)];

    for iGamma = 1:length(gamma)
        for iCost = 1:length(cost)
            s = ['-c ',num2str(cost(iCost)) ' -g ',num2str(gamma(iGamma)) ,' -b 1 -v 8']
            model(iGamma,iCost) = svmtrain(y, x, s) ;
        end

    end
end
toc


