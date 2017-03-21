% clear all
% dbstop if error
PLOT.boundary = 0;
PLOT.test = 0;
RLSC = 1
SVM = 1
mu{1} = [1 1];
mu{2} = [-1 -1];
S = [0.6 0.3; 0.3 0.6];
nTrain = 30;
nTest = 30;
colors = 'rgbk'
markers = 'oxox'
OUTLIER = round(nTrain*0.05)%10;

omu{1} = [-3 -1];
omu{2} = [3 1];
oS = [0.1 0.1; 0.1 0.1];

%%
figure(1),clf, hold on
for i = 1:2
    train{i} = mvnrnd(mu{i},S,nTrain);
    if OUTLIER
        train{i}(1:OUTLIER,:) = mvnrnd(omu{i},S,OUTLIER);
    end

    plot(train{i}(:,1),train{i}(:,2),'color',colors(i),'marker',markers(i),'linestyle','none')
end
for i = 1:2
    test{i} = mvnrnd(mu{i},S,nTest);
    if PLOT.test
        plot(test{i}(:,1),test{i}(:,2),'color',colors(i),'marker',markers(i),'linestyle','none')
    end
end
axis tight
vx = xlim
vy = ylim

%%
y = [ones(nTrain,1);zeros(nTrain,1)];
testY = [ones(nTest,1);zeros(nTest,1)];

clear corrTest corrTrain w
scale = 'log'
switch scale
    case 'linear'
        vLambda = linspace(0,100,3)
    case 'log'
        vLambda = [0 logspace(-3,1,4)]
end

varsToUse = [];

%%
for iLambda = 1:length(vLambda)
    lambda = vLambda(iLambda);
    x = [train{1};train{2}];
    testX = [test{1};test{2}];

    if RLSC
        w(iLambda,:) = trainRLSC(x,y,varsToUse,lambda);

        if PLOT.boundary
            colors = 'bk';
            markers = 'ox';
            % clear h
            h(1) = plot(vx,-vx*w(iLambda,2)/w(iLambda,3)-w(iLambda,1)/w(iLambda,3),'g-','linewidth',iLambda);
            h(2) = plot(vx,-vx*w(iLambda,2)/w(iLambda,3)-w(iLambda,1)/w(iLambda,3)+1/w(iLambda,3),'r-','linewidth',iLambda);
            h(3) = plot(vx,-vx*w(iLambda,2)/w(iLambda,3)-w(iLambda,1)/w(iLambda,3)+.5/w(iLambda,3),'b-','linewidth',iLambda);
            pause
        end
        %%
        % test error
        trainY = w(iLambda,:)*[ones(size(x,1),1), x]';
        corrTrain(iLambda) = ( sum(trainY(1:nTrain)>0.5) +  sum(trainY(nTrain+1:2*nTrain)<=0.5) )/nTrain/2;
        %%
        predY = w(iLambda,:)*[ones(size(testX,1),1), testX]';
        predY(1:nTest) = predY(1:nTest)>0.5;
        predY(nTest+1:2*nTest) = predY(nTest+1:2*nTest)<=0.5;
        
        corrTest(iLambda) = mean(predY);  % ( sum(testY(1:nTest)>0.5) +  sum(testY(nTest+1:2*nTest)<=0.5) ) /nTest/2 ;
        disp(['lambda =',num2str(lambda)])
        disp(['weight =' num2str(w(iLambda,:))])
        disp([num2str(corrTest(iLambda)) ,'/' ,num2str(corrTrain(iLambda)), '(test/train)'])
    end
    if SVM

        model = svmtrain(y, x, '-c 1 -g 2 -b 1');
        % model = svmtrain(heart_scale_label, heart_scale_inst, '-c 1 -g 2 -b 1')
        %        w(iLambda,:) = trainRLSC(x,y,varsToUse,lambda);

%                 if PLOT.boundary
%                     colors = 'bk';
%                     markers = 'ox';
%                     % clear h
%                     h(1) = plot(vx,-vx*w(iLambda,2)/w(iLambda,3)-w(iLambda,1)/w(iLambda,3),'g-','linewidth',iLambda);
%                     h(2) = plot(vx,-vx*w(iLambda,2)/w(iLambda,3)-w(iLambda,1)/w(iLambda,3)+1/w(iLambda,3),'r-','linewidth',iLambda);
%                     h(3) = plot(vx,-vx*w(iLambda,2)/w(iLambda,3)-w(iLambda,1)/w(iLambda,3)+.5/w(iLambda,3),'b-','linewidth',iLambda);
%                     pause
%                 end
        %%
        %%
        % test error
        %        trainY = w(iLambda,:)*[ones(size(x,1),1), x]';
        %        corrTrain(iLambda) = ( sum(trainY(1:nTrain)>0.5) +  sum(trainY(nTrain+1:2*nTrain)<=0.5) )/nTrain/2;
        %%
        %       testX = [test{1};test{2}];
        %        testY = w(iLambda,:)*[ones(size(testX,1),1), testX]';

        [tmp, corrTrainSVM(:,iLambda), decision_values] = svmpredict(y, x , model) ;
        [predY, corrTestSVM(:,iLambda), decision_values] = svmpredict(testY, testX , model) ;
        % [predicted_label, accuracy, decision_values/prob_estimates] = svmpredict(testing_label_vector, testing_instance_matrix, model [,'libsvm_options']);



        % corrTestSVM(iLambda) = mean(predY==testY); % ( sum(testY(1:nTest)>0.5) +  sum(testY(nTest+1:2*nTest)<=0.5) ) /nTest/2 ;
        disp(['lambda =',num2str(lambda)])
        disp(['weight =' num2str(w(iLambda,:))])
        disp([num2str(corrTestSVM(1,iLambda)) ,'/' ,num2str(corrTrain(iLambda)), '(test/train)'])
    end

end
%%
figure(2),clf,hold on,
vxLambda = vLambda
if vLambda(1)==0
    vxLambda(1) = vLambda(2)^2/vLambda(3)
end
plot(vxLambda,corrTrain,'bx:')
plot(vxLambda,corrTest,'bo-')
plot(vxLambda,corrTrainSVM(1,:)/100,'rx:')
plot(vxLambda,corrTestSVM(1,:)/100,'ro-')

legend('RLSC train','RLSC test','SVM train','SVM test')

set(gca,'xscale',scale)
