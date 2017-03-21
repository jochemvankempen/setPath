function [P] = decodeRLSC(cResp,DEC)

%%

if size(cResp,1) == 2
        [ind_train,ind_test]    = leaveOneOutRLSC_mem(cResp,DEC);
        [P]                     = evaluatePerformance_mem(cResp,DEC,ind_train,ind_test);
                
elseif size(cResp,1) >= 2
    %         error('needs fixing')
    %         [P,CC] = leaveOneOutRLSC_OVA(cResp,DEC);
    %         [ind_train,ind_test]    = leaveOneOutRLSC_mem_OVA(cResp,DEC);
    error
    [ind_train,ind_test]    = leaveOneOutRLSC_mem(cResp,DEC);
    [P]                     = evaluatePerformance_mem_OVA(cResp,DEC,ind_train,ind_test);
    %         [P,CC] = evaluatePerformance(P,CC);
    
end
%     P(iTime,:,:)=x;
