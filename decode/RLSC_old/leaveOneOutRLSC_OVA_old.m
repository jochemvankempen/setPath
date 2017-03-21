function [P,C] = leaveOneOutRLSC_OVA(cResp,C)
% WARNING INCOMPLETE - OBSOLETE - TO BE DOUBLE-CHECKED
%% 07 Feb 18
% one versus all classifier
%% 07 Feb 10
% Use RLSC for classification
% Compute performance for a range of lambda values.
% Decide which lambda at the end, by looking at the best test performance
% use all voxels

%% get the number of trials per session, index for each trial
nSession = length(C.vSession);

%% compute the range of lambda

allResp = [];
for iClass = 1:length(cResp)
    allResp = [allResp;cResp{iClass}];
end

%%
tmp = allResp*allResp';
cLambda = round(log10(max(tmp(:))));
vLamRange = [-5 2];
nLam = 5 ; % 15;
vLambda = logspace( cLambda + vLamRange(1) , cLambda + + vLamRange(2), nLam);
C.cLambda = cLambda;
C.vLamRange = vLamRange;
C.nLam = nLam;
C.vLambda = vLambda;

for iSample = 1:nSession
    if mod(iSample,20)==0
        disp([num2str(iSample),'th sample'])
    end

    % select a leave-one-out test and training
    vTrial = [1:12] +(iSample-1)*12;
    for iClass = 1:length(cResp)
        test{iClass}= cResp{iClass}(vTrial,:);
        train{iClass} = cResp{iClass}(setdiff([1:size(cResp{iClass},1)],vTrial),:);
    end

    %% train and test 4 classifiers with one-versus-all scheme
    [nClTrial(1)] = size(train{1},1);
    [nClTrial(2)] = size(train{2},1);
    [nTsTrial(1)] = size(test{1},1);
    [nTsTrial(2)] = size(test{2},1);

    x = []; testX = [];
    for iCl = 1:length(cResp)
        x = [x;train{iCl}];
        testX = [testX;test{iCl};];
    end

    for iClassifier = 1:length(cResp)
        y = [];
        
        for iCl = 1:length(cResp)
            if iCl == iClassifier
                y = [y;ones(nClTrial(1),1);];
            else
                y = [y;zeros(nClTrial(1),1);];
            end
        end


        for iLambda = 1:length(vLambda)
            lambda = vLambda(iLambda);
            w = trainRLSC(x,y,[],lambda);
            trainY = w*[ones(size(x,1),1),x]';
            predY = w*[ones(size(testX,1),1), testX]';

            P(iClassifier,iLambda,iSample).w = w;
            P(iClassifier,iLambda,iSample).trainY = trainY;
            P(iClassifier,iLambda,iSample).predY = predY;
        end
    end

    

    %         cl = cell(1,2);
    %         ts = cell(1,2);
    %         for iCl = 1:4
    %             if iCl == iClassifier
    %                 cl{1} = train{iCl};
    %                 ts{1} = test{iCl};
    %             else
    %                 cl{2} = [cl{2};train{iCl}];
    %                 ts{2} = [ts{2};test{iCl}];
    %             end
    %         end
    %         [nClTrial(1)] = size(cl{1},1);
    %         [nClTrial(2)] = size(cl{2},1);
    %         [nTsTrial(1)] = size(ts{1},1);
    %         [nTsTrial(2)] = size(ts{2},1);
    %
    %         C.nClTrial = nClTrial;
    %         C.nTsTrial = nTsTrial;
    %
    %         y = [ones(nClTrial(1),1);zeros(nClTrial(2),1)];
    %         x = [cl{1};cl{2}];
    %         testX = [ts{1};ts{2}];
    %
    %         for iLambda = 1:length(vLambda)
    %             lambda = vLambda(iLambda);
    %             w = trainRLSC(x,y,[],lambda);
    %             trainY = w*[ones(size(x,1),1),x]';
    %             predY = w*[ones(size(testX,1),1), testX]';
    %
    %             P(iClassifier,iLambda,iSample).w = w;
    %             P(iClassifier,iLambda,iSample).trainY = trainY;
    %             P(iClassifier,iLambda,iSample).predY = predY;
    %         end
    %   end
end

C.nClTrial = nClTrial;
C.nTsTrial = nTsTrial;
