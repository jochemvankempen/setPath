nSimulation = 100;
%% 

allRLSC = [];
allSVM = [];

for iSimulation = 1:nSimulation
    demoRLSC 
    allRLSC(iSimulation,:) = corrTest;
    allSVM(iSimulation,:) = corrTestSVM(1,:)/100;
end
%% 
xshift = 1.2
figure(3),cla, hold on 
plot(vxLambda./xshift,mean(allRLSC),'bo-')
plot([vxLambda;vxLambda]./xshift,prctile(allRLSC,[2.5 97.5]),'b-')
plot(vxLambda.*xshift,mean(allSVM),'ro-')
plot([vxLambda;vxLambda].*xshift,prctile(allSVM,[2.5 97.5]),'r-')


legend('RLSC test','SVM test')
set(gca,'xscale',scale)
