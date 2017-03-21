%% leave one out validation
% svm
% strange why the performance doesn't drop to .5 for 1 voxel
% n-fold validation 
% create test and training at once 

clear all 
dbstop if error 
PLOT.covariance = 0;
PLOT.ttest = 0;

nVoxel = 100 % 1000;
separation = .1;
mu{1} = separation*randn(nVoxel,1).*ones(nVoxel,1);
mu{2} = -separation*randn(nVoxel,1).*ones(nVoxel,1);
correlation = 0.1;
% S = correlation*ones(nVoxel)+eye(nVoxel)*(1-correlation);
 S = 0.05*rand(nVoxel)+eye(nVoxel);
 S = S'*S;

if PLOT.covariance
    figure(3),
    imagesc(S),colorbar
end

nTrain = 50;
nTest = 50;
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
    train{i} = mvnrnd(mu{i},S,nTrain);
    if OUTLIER
        train{i}(1:OUTLIER,:) = mvnrnd(omu{i},S,OUTLIER);
    end

    plot(train{i}(:,1),train{i}(:,2),'color',colors(i),'marker',markers(i),'linestyle','none')
end
colors = 'bk';
markers = 'ox';
PLOT.test = 1;

for i = 1:2
    test{i} = mvnrnd(mu{i},S,nTest);
    if OUTLIER
        test{i}(1:OUTLIER,:) = mvnrnd(omu{i},S,OUTLIER);
    end
    if PLOT.test
        plot(test{i}(:,1),test{i}(:,2),'color',colors(i),'marker',markers(i),'linestyle','none')
    end
end
axis tight
vx = xlim
vy = ylim

%%
y = [ones(nTrain,1);zeros(nTrain,1)];

PLOT.boundary = 0;
saveW = 0 ;
clear corrTest corrTrain w
scale = 'log'
nLambda = 3
switch scale
    case 'linear'
        vLambda = linspace(0,1,nLambda)
    case 'log'
        vLambda = [0 logspace(-nLambda+1,0,nLambda)]
end

% keyboard
% vUsedVoxel = 1:10:100% nVoxel
vUsedVoxel = 1:10:nVoxel
% vUsedVoxel = [1:10:100 200 400 1000] %
% 1:200:nVoxel
for iLambda = 1:length(vLambda)
    lambda = vLambda(iLambda)
    for iVoxel = vUsedVoxel
        % [h,p,ci,stat] = ttest(train{1},train{2});
        % [sortTstat,iSortTstat] = sort(abs(stat.tstat));
        % varsToUse = iSortTstat(1:iVoxel);
        varsToUse = randperm(nVoxel);
        varsToUse = varsToUse(1:iVoxel);
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
        trainY = w*[ones(size(x,1),1), x(:,varsToUse)]';
        corrTrain(iLambda,iVoxel) = ( sum(trainY(1:nTrain)>0.5) +  sum(trainY(nTrain+1:2*nTrain)<=0.5) )/nTrain/2;
        %    disp(corrTrain/nTrain/2)
        %%
        testX = [test{1};test{2}];
        testY = w*[ones(size(testX,1),1), testX(:,varsToUse)]';
        corrTest(iLambda,iVoxel) = ( sum(testY(1:nTest)>0.5) +  sum(testY(nTest+1:2*nTest)<=0.5) ) /nTest/2 ;
        %    disp(corrTest/nTest/2)
        disp([num2str(iVoxel),' voxels: ' ,num2str(corrTest(iLambda,iVoxel)) ,'/' ,num2str(corrTrain(iLambda,iVoxel)), '(test/train)'])
        if saveW
            W{iVoxel}(iLambda,:) = w;
        end
    end
end
%%
figure(2),clf,hold on,
colors = 'rgbkymc'
for iLambda = 1:length(vLambda)
    plot(vUsedVoxel,corrTrain(iLambda,vUsedVoxel),'x:','color',colors(iLambda))
    legendtext{iLambda} = ['lm=',num2str(vLambda(iLambda))];
end
legend(legendtext,'location','northwest')
for iLambda = 1:length(vLambda)
    plot(vUsedVoxel,corrTest(iLambda,vUsedVoxel),'o-','color',colors(iLambda))
end
%plot(vLambda,corrTrain,'bo-')
%plot(vLambda,corrTest,'kx-')
%%
set(gca,'xscale',scale)
%%
if PLOT.ttest
    [h,p,ci,stat] = ttest(train{1},train{2});
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