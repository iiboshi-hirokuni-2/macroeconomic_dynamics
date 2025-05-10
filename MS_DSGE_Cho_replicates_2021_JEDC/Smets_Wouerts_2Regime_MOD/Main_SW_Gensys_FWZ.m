clc

addpath('./TZMatlabPrograms/')
addpath('./Func/')

clear all
 close all

 disp( '  '  )
disp('============================================')
disp('           start  Smets & Wouters (AER 2007) DSGE model ')
disp(' Solution method ;  Gensys version ')
disp('============================================')

 SW_model_gensys
 
%=====================
common_params_SW;

%% regime 1 (PMP/AFP) 
disp('------------------------------')
disp( 'regime 1 (PMP/AFP) '   )
disp('------------------------------')
% set parameter values to matrix 

Regime_parameters_R1;

% chek if all f evaluated at ss value are zero
check = eval(System_of_Eq)

GAM01 = eval(GAM0j);
GAM11 = eval(GAM1j);
PSI01 = eval(PSI0j);
PPI1 = eval(PPIj);

%  solve a DSGE model using GENSYS
[T1,TC,T0,fmat,fwt,ywt,gev,RC] = gensys(GAM01,GAM11,CC,PSI01,PPI1,1);
RC

T0_R1=T0; T1_R1 = T1;

%% regime 2 (AMP/PFP)
 h =2;
disp('------------------------------')
disp( 'regime 2 (AMP/PFP) '   )
disp('------------------------------')

% p_st(:,1)= p_st(:,i);
Regime_parameters_R2;

% chek if all f evaluated at ss value are zero
% check = eval(System_of_Eq)

GAM02 = eval(GAM0j);
GAM12 = eval(GAM1j);
PSI02 = eval(PSI0j);
PPI2  = eval(PPIj);

[T1,TC,T0,fmat,fwt,ywt,gev,RC] = gensys(GAM02,GAM12,CC,PSI02,PPI2,1);
RC

T0_R2=T0; T1_R2 = T1;


%%   FWZ method for BM AER model
  SW_FWZ_methods_R2
 