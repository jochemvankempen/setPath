%% demo linear projection
% based on papers 10.1038/78856, 10.1073/pnas.1317557111 and 10.1111/ejn.12859
% jochem van kempen 22/02/2017


figure(1),clf
vector = 1:0.01:3; %make one template vector

nTrials=100;
subplot(2,2,1)
plot(vector)
title('single vector')

vectors = zeros(length(vector),nTrials);
for ivec = 1:nTrials
    %make several vectors with different slopes and add noise
    vectors(:,ivec) = vector*rand + rand(1,length(vector))/10;
end

subplot(2,2,2)
hold on
plot(vectors)
plot(mean(vectors,2),'k','linewidth',5)
title('all vectors with added noise')


averageVector   = mean(vectors,2); %column vector, average across trials
normVector      = averageVector/sqrt(averageVector'*averageVector)^2;% norm of the average vector
subplot(2,2,3)
plot(normVector)
title('norm')

lp = zeros(nTrials,1);
for itrial = 1:nTrials
    % for each trial, compute linear projection by multiplying single trial
    % row vector by the norm and deviding by the length
    lp(itrial,1) = (vectors(:,itrial)'*normVector)/length(normVector);
    
end
subplot(2,2,4)
plot(lp, vectors(end,:)','.')
title('relation endpoint vector, projection')



