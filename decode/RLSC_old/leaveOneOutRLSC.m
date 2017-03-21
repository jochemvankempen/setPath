function [P,C] = leaveOneOutRLSC(cResp,C,trialnum)
% 07 Mar 17
% changed the class assignment strategy -- change the # of test and
% training according to the # of maximally available trials
%
% nClTrial( >> nClTrial(iSample,
% nTsTrial( >> nTsTrial(iSample,

%% 07 Feb 10
% Use RLSC for classification
% Compute performance for a range of lambda values.
% Decide which lambda at the end, by looking at the best test performance
% use all voxels

%% get the number of sessions -- same as n-Fold cross validation by default
% nSession = length(C.vSession);

%% compute the range of lambda

allResp = [];
for iClass = 1:length(cResp)
    allResp = [allResp;cResp{iClass}];
end

%%
if ~isfield(C,'nLam')
    C.nLam = 5 ;
end
if ~isfield(C,'vLamRange')
    C.vLamRange = [-5 2];
end
if ~isfield(C,'vLambda')
    tmp = allResp*allResp';
    C.cLambda = round(log10(max(tmp(:))));
    if C.nLam == 1
        C.vLambda = logspace( C.cLambda + C.vLamRange(1) , C.cLambda + + C.vLamRange(1), C.nLam);
    else
        C.vLambda = logspace( C.cLambda + C.vLamRange(1) , C.cLambda + + C.vLamRange(2), C.nLam);
    end
else
    C.nLam = length(C.vLambda);
    % C.vLamRange = nan;
end

if ~isfield(C,'nFoldValidation')
    C.nFoldValidation = C.vSession;
end
if ~isfield(C,'randomSampleTest')
    C.randomSampleTest = []; % Gabriel used 30% for test, 70% for training.
    % this may be necessary to get better performance when the # of trials
    % are uneven (like Kanai's data, Alan's data, Hackjin's data)
end
%%
% select a leave-one-out test and training
% create test and trial trials for nFoldValidation
if 0
    % old method
    % specific to face/house task
    vTrial = [1:12] +(iSample-1)*12;
    for iClass = 1:4
        test{iClass}= cResp{iClass}(vTrial,:);
        train{iClass} = cResp{iClass}(setdiff([1:size(cResp{iClass},1)],vTrial),:);
    end
end
if 0
    % 07 Mar 14
    % new method -- works for unbalanced data -- but multiple of
    % nFoldValidation
    for iClass = 1:length(cResp)
        nTrialPerSession(iClass) = size(cResp{iClass},1);
        nTrialPerValidation(iClass) = nTrialPerSession(iClass)/C.nFoldValidation;
    end
    for iClass = 1:length(cResp)
        vTestTrial = [1:nTrialPerValidation(iClass)] +(iSample-1)*nTrialPerValidation(iClass);
        test{iClass}= cResp{iClass}(vTestTrial,:);
        vTrainTrial = setdiff([1:size(cResp{iClass},1)],vTestTrial);
        train{iClass} = cResp{iClass}(vTrainTrial,:);
    end
end


if isempty(C.randomSampleTest)

    % 07 Mar 17 new method -- works for arbitrary # of trials for each class
    % count the number of trials for each validation set
    for iClass = 1:length(cResp)
        nTrialPerSession(iClass) = size(cResp{iClass},1);
        nCeiledTrialPerValidation(iClass) = ceil(nTrialPerSession(iClass)/C.nFoldValidation);
        nFlooredSession(iClass) = nCeiledTrialPerValidation(iClass)*C.nFoldValidation -  nTrialPerSession(iClass) ;
    end

    jEndLastTestTrial = zeros(length(cResp));
    for iSample = 1:C.nFoldValidation
        for iClass = 1:length(cResp)
            if iSample <= nFlooredSession(iClass)
                nTrialPerValidation(iClass,iSample) = floor(nTrialPerSession(iClass)/C.nFoldValidation);
            else
                nTrialPerValidation(iClass,iSample) = ceil(nTrialPerSession(iClass)/C.nFoldValidation);
            end
        end

        for iClass = 1:length(cResp)
            vTestTrial{iSample,iClass} = [1:nTrialPerValidation(iClass,iSample)] + jEndLastTestTrial(iClass);
            jEndLastTestTrial(iClass) =  vTestTrial{iSample,iClass}(end);
            test{iSample,iClass}= cResp{iClass}(vTestTrial{iSample,iClass},:);
            vTrainTrial{iSample,iClass} = setdiff([1:nTrialPerSession(iClass)],vTestTrial{iSample,iClass});
            train{iSample,iClass} = cResp{iClass}(vTrainTrial{iSample,iClass},:);

            if length(unique([vTestTrial{iSample,iClass},vTrainTrial{iSample,iClass}]) ) ~= nTrialPerSession (iClass)
                error('# of test and training trials doesn''t match # of trials for this session')
            end
        end
    end
else
    if C.randomSampleTest == -1; % leave one pair out 
        keyboard 
    elseif C.randomSampleTest == 0
        % 08 Mar 26
        % to cope with drifting neuronal response problem, drift the test
        % and training accordingly 
 
        if 1
            % keyboard
            % find the closest N trial for test and training 
            if length(cResp) ~= 2 
                error('')
            end
            C.nFoldValidation = size(cResp{1},1) + size(cResp{2},1);
            
            % create sample for left trial 
            jSample = 0; 
            for iClass = 1:length(cResp)
                
                for iSample = 1:size(cResp{iClass},1)
                    jSample = jSample + 1; 
                    
                    test{jSample,3-iClass} = []; % no test for the other class 
                    test{jSample,iClass}   = cResp{iClass}(iSample,:);
                    
                    % find the closest trials for training 
                    testTrialNum(jSample) = trialnum{iClass}(iSample);
                    
                    % own class 
                    [tmp,iOwnTmp{jSample}] = sort( abs ( trialnum{iClass} - testTrialNum(jSample) ));
                    vOwnTmp{jSample} = sort(iOwnTmp{jSample}(2:C.nHalfTraining*2+1));
                    vOwnTrialnum(jSample,:) = trialnum{iClass}(vOwnTmp{jSample});
                    train{jSample,iClass} = cResp{iClass}(vOwnTmp{jSample} ,:);
                    % the other class 
                    [tmp,iOtherTmp{jSample}] = sort( abs ( trialnum{3-iClass} - testTrialNum(jSample) ));
                    vOtherTmp{jSample} = sort(iOtherTmp{jSample}(1:C.nHalfTraining*2));;
                    vOtherTrialnum(jSample,:) = trialnum{3-iClass}(vOtherTmp{jSample});
                    train{jSample,3-iClass} = cResp{3-iClass}(vOtherTmp{jSample},:);
                    
                end
            end
                if 0  % sanity check 
                    testTrialNum
                vOwnTrialnum
                vOtherTrialnum
                end 
            
            

        else%% this method still has a drifting problem, especially towards end.  Use trialnum to avoid this problem
            % 1. sample random number of equale trials (for test and training)
            for iClass = 1:length(cResp)
                nTrialPerSession(iClass) = size(cResp{iClass},1);
            end
            nTest = 1 ;
            nTrain = C.nHalfTraining * 2;

            C.nFoldValidation = min(nTrialPerSession);

            if nTrain + nTest > C.nFoldValidation
                disp(' # of 1 + 2*C.nHalfTraining > min(nTrialPerSession) ')
            end

            for iSample = 1:C.nFoldValidation
                for iClass = 1:length(cResp)

                    vTestTrial{iSample,iClass} = iSample;
                    if iSample - 1 < C.nHalfTraining
                        vTrainTrial{iSample,iClass} = [1:(iSample-1) , (iSample+1):(C.nHalfTraining*2+1)] ;
                    elseif iSample + C.nHalfTraining > C.nFoldValidation
                        vTrainTrial{iSample,iClass} = [(C.nFoldValidation - 2*C.nHalfTraining ):(iSample-1) ,...
                            (iSample+1):(C.nFoldValidation)] ;
                    else
                        vTrainTrial{iSample,iClass} = [(iSample-C.nHalfTraining):(iSample-1) ,...
                            (iSample+1):(iSample+C.nHalfTraining)] ;
                    end
                    test {iSample,iClass} = cResp{iClass}(vTestTrial {iSample,iClass},:);
                    train{iSample,iClass} = cResp{iClass}(vTrainTrial{iSample,iClass},:);
                end
            end
        end
        % keyboard

    else

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
    end
end

%%

% for iSample = 1:nSession
for iSample = 1:C.nFoldValidation
    if mod(iSample,20)==0 & 0 
        disp([num2str(iSample),'th fold validation'])
    end

    cl = [];
    switch C.task
        case {'FaceCheck','HappyFear','TheoryOfMind','AVsync','Varela',...
                'sfm','Gender','SacNosac','UpInv','FirstSecond','AttendNonattend','erp',...
                'fastSlowRT_CFS','CFSvsFS','cfsbm',...
                'IowaLocalizer_faceVsAll','IowaLocalizer_faceVsMondrian',...
                'TW_IowaLocalizer','surrogate','TW_IowaCFS','PAC_IowaCFS'}
            cl{1} = train{iSample,1};
            cl{2} = train{iSample,2};
            ts{1} = test{iSample,1};
            ts{2} = test{iSample,2};
        case 'FaceHouse'
            cl{1} = [train{iSample,1};train{iSample,2}];
            cl{2} = [train{iSample,3};train{iSample,4}];
            ts{1} = [test{iSample,1};test{iSample,2}];
            ts{2} = [test{iSample,3};test{iSample,4}];
        case 'FearNeut'
            cl{1} = [train{iSample,1};train{iSample,3}];
            cl{2} = [train{iSample,2};train{iSample,4}];
            ts{1} = [test{iSample,1};test{iSample,3}];
            ts{2} = [test{iSample,2};test{iSample,4}];
        case 'AttendFearNeut'
            cl{1} = [train{iSample,1}];
            cl{2} = [train{iSample,2}];
            ts{1} = [test{iSample,1}];
            ts{2} = [test{iSample,2}];
        case 'IgnoreFearNeut'
            cl{1} = [train{iSample,3}];
            cl{2} = [train{iSample,4}];
            ts{1} = [test{iSample,3}];
            ts{2} = [test{iSample,4}];
    end

    [nClTrial(iSample,1)] = size(cl{1},1);
    [nClTrial(iSample,2)] = size(cl{2},1);
    [nTsTrial(iSample,1)] = size(ts{1},1);
    [nTsTrial(iSample,2)] = size(ts{2},1);

    C.nClTrial = nClTrial;
    C.nTsTrial = nTsTrial;

    y = [ones(nClTrial(iSample,1),1);-ones(nClTrial(iSample,2),1)];
    %    y = [ones(nClTrial(iSample,1),1);zeros(nClTrial(iSample,2),1)];
    x = [cl{1};cl{2}];
    testX = [ts{1};ts{2}];

    switch C.spaceAverage
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
    for iLambda = 1:length(C.vLambda)
        lambda = C.vLambda(iLambda);
        w = trainRLSC(x,y,[],lambda);
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

        P(iLambda,iSample).w = w;
        P(iLambda,iSample).trainY = trainY;
        P(iLambda,iSample).predY = predY;
        % P(iLambda,iSample).corrTrain = corrTrain;
        % P(iLambda,iSample).corrTest = corrTest;

        %         if corrTrain <= corrTest & corrTrain < 1
        %             disp('too much penalty')
        %             break
        %         end


    end
end


