
function [ RC_1, RC_2]=func_Model_SW_gsys(par)

SW_model_gensys

%=====================
common_params_SW;

%% regime 1 (PMP/AFP) 
% nst =1; i=1;
Regime_parameters_R1;
crpi = par(1);

% chek if all f evaluated at ss value are zero
% check = eval(System_of_Eq)

GAM01 = eval(GAM0j);
GAM11 = eval(GAM1j);
PSI01 = eval(PSI0j);
PPI1 = eval(PPIj);

%  solve a DSGE model using GENSYS
[T1,TC,T0,fmat,fwt,ywt,gev,RC_1] = gensys(GAM01,GAM11,CC,PSI01,PPI1,1);
RC_1 = RC_1';

%% regime 2 (AMP/PFP)
% i=2;
% p_st(:,1)= p_st(:,i);
Regime_parameters_R2;
crpi = par(2);

% chek if all f evaluated at ss value are zero
% check = eval(System_of_Eq)

GAM02 = eval(GAM0j);
GAM12 = eval(GAM1j);
PSI02 = eval(PSI0j);
PPI2  = eval(PPIj);

[T1,TC,T0,fmat,fwt,ywt,gev,RC_2] = gensys(GAM02,GAM12,CC,PSI02,PPI2,1);
RC_2 = RC_2';

RC=(RC_1.*RC_2)';
 








