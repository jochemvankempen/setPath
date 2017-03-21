function [P] = decode_SVM(data,DEC)

[ind_train,ind_test]    = decode_permTrialIdx(data,DEC);

accuracy = zeros(1,DEC.nCrossVal);
for iSample = 1:DEC.nCrossVal
    if mod(iSample,10)==0
        disp([num2str(iSample),'th fold validation'])
    end
    
    tr = cell(1,length(data));
    ts = cell(1,length(data));
    nTrTrial = zeros(1,length(data));
    nTsTrial = zeros(1,length(data));
    for iCat = 1:length(data)
        tr{iCat} = data{iCat}(squeeze(ind_train(iSample,iCat,:)),:);
        ts{iCat} = data{iCat}(squeeze(ind_test(iSample,iCat,:)),:);
        [nTrTrial(iCat)] = size(tr{iCat},1);
        [nTsTrial(iCat)] = size(ts{iCat},1);
    end
    [tr,ts]=normalizePatternsEachFold(tr,ts);
    
    DEC.nTrTrial = nTrTrial;
    DEC.nTsTrial = nTsTrial;
    
    yTrain  = [];
    yTest   = [];
    xTrain  = [];
    xTest   = [];
    for iCat = 1:length(data)
        
        yTrain = [yTrain ones(nTrTrial(iCat),1)*iCat];
        yTest  = [yTest ones(nTsTrial(iCat),1)*iCat];
        
        xTrain = cat(1,xTrain,tr{iCat});
        xTest  = cat(1,xTest,ts{iCat});
    end
    
    
    %%
    
    if length(data) == 2
        yTrain  = abs(yTrain-2);
        yTest   = abs(yTest-2);
        
        [svmModel]             = svmtrain(yTrain, xTrain, DEC.trnParams);
        [~, tmpAccuracy, ~]    = svmpredict(yTest, xTest, svmModel, DEC.tstParams);
    elseif length(data) > 2
        svmModel               = ovrtrain(yTrain, xTrain, DEC.trnParams);
        [~, tmpAccuracy, ~]    = ovrpredict(yTest, xTest, svmModel, DEC.tstParams);
    end
        
    accuracy(iSample) = tmpAccuracy(1);
end

% P.predLabel         = predict_label;
P.testAccuracy      = mean(accuracy);
P.testAccuracyStd   = std(accuracy);

%train the SVM with options as specified
 
% `svm-train' Usage
% =================
% matlab> model = svmtrain(training_label_vector, training_instance_matrix [, 'libsvm_options']);
% 
%         -training_label_vector:
%             An m by 1 vector of training labels (type must be double).
%         -training_instance_matrix:
%             An m by n matrix of m training instances with n features.
%             It can be dense or sparse (type must be double).
%         -libsvm_options:
%             A string of training options in the same format as that of LIBSVM
% 
% options:
% -s svm_type : set type of SVM (default 0)
% 	0 -- C-SVC		(multi-class classification)
% 	1 -- nu-SVC		(multi-class classification)
% 	2 -- one-class SVM	
% 	3 -- epsilon-SVR	(regression)
% 	4 -- nu-SVR		(regression)
% -t kernel_type : set type of kernel function (default 2)
% 	0 -- linear: u'*v
% 	1 -- polynomial: (gamma*u'*v + coef0)^degree
% 	2 -- radial basis function: exp(-gamma*|u-v|^2)
% 	3 -- sigmoid: tanh(gamma*u'*v + coef0)
% 	4 -- precomputed kernel (kernel values in training_set_file)
% -d degree : set degree in kernel function (default 3)
% -g gamma : set gamma in kernel function (default 1/num_features)
% -r coef0 : set coef0 in kernel function (default 0)
% -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
% -n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
% -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
% -m cachesize : set cache memory size in MB (default 100)
% -e epsilon : set tolerance of termination criterion (default 0.001)
% -h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)
% -b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
% -wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)
% -v n: n-fold cross validation mode
% -q : quiet mode (no outputs)
% 
% 
% The k in the -g option means the number of attributes in the input data.
% 
% option -v randomly splits the data into n parts and calculates cross
% validation accuracy/mean squared error on them.
% 
% See libsvm FAQ for the meaning of outputs.
% 
% Returned Model Structure
% ========================
% 
% The 'svmtrain' function returns a model which can be used for future
% prediction.  It is a structure and is organized as [Parameters, nr_class,
% totalSV, rho, Label, ProbA, ProbB, nSV, sv_coef, SVs]:
% 
%         -Parameters: parameters
%         -nr_class: number of classes; = 2 for regression/one-class svm
%   [predicted_label, accuracy, prob_estimates] = libsvmpredict(tstLabels, tstData, SVM_Model,  libSVM_TstParams);      -totalSV: total #SV
%         -rho: -b of the decision function(s) wx+b
%         -Label: label of each class; empty for regression/one-class SVM
%         -ProbA: pairwise probability information; empty if -b 0 or in one-class SVM
%         -ProbB: pairwise probability information; empty if -b 0 or in one-class SVM
%         -nSV: number of SVs for each class; empty for regression/one-class SVM
%         -sv_coef: coefficients for SVs in decision functions
%         -SVs: support vectors
% 
% If you do not use the option '-b 1', ProbA and ProbB are empty
% matrices. If the '-v' option is specified, cross validation is
% conducted and the returned model is just a scalar: cross-validation
% accuracy for classification and mean-squared error for regression.
% 
% More details about this model can be found in LIBSVM FAQ
% (http://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html) and LIBSVM
% implementation document
% (http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf).


% 
% Usage: 
% [predicted_label, accuracy, decision_values/prob_estimates] = svmpredict(testing_label_vector, testing_instance_matrix, model [, 'libsvm_options']);

% options:
% -b probability_estimates: whether to predict probability estimates, 0 or 1 (default 0); for one-class SVM only 0 is supported
 
 
% Result of Prediction
% ====================
% 
% The function 'svmpredict' has three outputs. The first one,
% predictd_label, is a vector of predicted labels. The second output,
% accuracy, is a vector including accuracy (for classification), mean
% squared error, and squared correlation coefficient (for regression).
% The third is a matrix containing decision values or probability
% estimates (if '-b 1' is specified). If k is the number of classes
% in training data, for decision values, each row includes results of 
% predicting k(k-1)/2 binary-class SVMs. For classification, k = 1 is a
% special case. Decision value +1 is returned for each testing instance,
% instead of an empty vector. For probabilities, each row contains k values
% indicating the probability that the testing instance is in each class.
% Note that the order of classes here is the same as 'Label' field
% in the model structure.

%cross validataion example

% [trainY trainX] = libsvmread('./dna.scale');
% [testY testX] = libsvmread('./dna.scale.t');
% model = ovrtrain(trainY, trainX, '-c 8 -g 4');
% [pred ac decv] = ovrpredict(testY, testX, model);
% fprintf('Accuracy = %g%%\n', ac * 100);
% 
% Conduct CV on a grid of parameters 
% bestcv = 0;
% for log2c = -1:2:3,
%   for log2g = -4:2:1,
%     cmd = ['-q -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
%     cv = get_cv_ac(trainY, trainX, cmd, 3);
%     if (cv >= bestcv),
%       bestcv = cv; bestc = 2^log2c; bestg = 2^log2g;
%     end
%     fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', log2c, log2g, cv, bestc, bestg, bestcv);
%   end
% end

%one vs rest example




