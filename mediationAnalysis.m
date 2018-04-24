function mediationAnalysis(IV, MED, DV, varargin)


%%% comments:
%%% MacKinnon, D., Fairchild, A., & Fritz, M. (2007). Mediation Analysis. Annual Review of Psychology, 58(Hebb 1966), 593–602. https://doi.org/10.1146/annurev.psych.58.110405.085542.Mediation
%%%
%%% Four steps are involved in the Baron and Kenny approach to establishing
%%% mediation. 
%%% 1. (criteria1): a significant relation of the independent variable to the dependent 
%%% variable is required in Equation 1. This is a normal regression, not
%%% taking into account any other variables (Y = 1 + cX + e1).
%%% 2. (criteria2): a significant relation of the independent variable to
%%% the hypothesized mediating variable is required in Equation 3. This is
%%% also a normal regression. (Y = 1 + aX + e3).
%%% 3. (criteria3). the mediating variable must be significantly related to
%%% the dependent variable when both the independent variable and mediating
%%% variable are predictors of the dependent variable in Equation 2. This
%%% is a dependent regression (Y = 1 + c'X + bM +e2).


%%% figure 1. 
disp('      M ')
disp('    /   \')
disp('  a/     \b')
disp('  /       \')
disp(' /         \')
disp('IV -- c'' -- DV')





