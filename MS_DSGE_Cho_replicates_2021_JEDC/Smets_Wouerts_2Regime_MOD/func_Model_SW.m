
function [ PP, A,B,C]=func_Model_SW(par)

  
 %% Load  Smets & Wouters (AER 2007) DSGE model
SW_model;

%%

B1j_R1 = -1*jacobian(System_of_Eq, Y_t);
A1j_R1 = jacobian(System_of_Eq, Y_f);
B2j_R1  = jacobian(System_of_Eq, Y_b);
C1j_R1 = jacobian(System_of_Eq, Epsilon_t);

B1j_R2 = -1*jacobian(System_of_Eq, Y_t);
A1j_R2 = jacobian(System_of_Eq, Y_f);
B2j_R2  = jacobian(System_of_Eq, Y_b);
C1j_R2 = jacobian(System_of_Eq, Epsilon_t);

%  transtion Probs
  MS_transition_prob;
  
  %=====================
common_params_SW;

% regime 1 (PMP/AFP) 
Regime_parameters_R1;

%========================
crpi     =  par(1);
cry      =    par(3);
%========================

B1{1} = eval(B1j_R1);
A1{1} = eval(A1j_R1);
B2{1} = eval(B2j_R1);
C1{1} = eval(C1j_R1);

%                 A=B1^(-1)*A1
%                 B=B1^(-1)*B2
%                 C=B1^(-1)*C1
i=1;
            A{i,1}=B1{i}\A1{i};
            B{i,1}=B1{i}\B2{i};
            C{i,1}=B1{i}\C1{i};  

% regime 2 (AMP/PFP) 
Regime_parameters_R2;

%========================
crpi      =  par(2);
cry        =  par(4);
%========================

B1{2} = eval(B1j_R1);
A1{2} = eval(A1j_R1);
B2{2} = eval(B2j_R1);
C1{2} = eval(C1j_R1);

%                 A=B1^(-1)*A1
%                 B=B1^(-1)*B2
%                 C=B1^(-1)*C1
i=2;
            A{i,1}=B1{i}\A1{i};
            B{i,1}=B1{i}\B2{i};
            C{i,1}=B1{i}\C1{i};  
 








