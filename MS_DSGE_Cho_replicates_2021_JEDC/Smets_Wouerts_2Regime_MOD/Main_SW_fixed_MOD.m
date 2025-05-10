% addpath('./TZMatlabPrograms/')
addpath(genpath('./MODMethod'))
addpath('./Func/')

clear all
close all

SW_model;

B1j_R1 = -1*jacobian(System_of_Eq, Y_t);
A1j_R1 = jacobian(System_of_Eq, Y_f);
B2j_R1  = jacobian(System_of_Eq, Y_b);
C1j_R1 = jacobian(System_of_Eq, Epsilon_t);

B1j_R2 = -1*jacobian(System_of_Eq, Y_t);
A1j_R2 = jacobian(System_of_Eq, Y_f);
B2j_R2  = jacobian(System_of_Eq, Y_b);
C1j_R2 = jacobian(System_of_Eq, Epsilon_t);

%=====================
common_params_SW;

%% regime 1 (PMP/AFP) 
disp('------------------------------')
disp( 'regime 1 (PMP/AFP) '   )
disp('------------------------------')
% set parameter values to matrix 

% p_st(:,1)= p_st(:,i);
 Regime_parameters_R1;

B1{1} = eval(B1j_R1);
A1{1} = eval(A1j_R1);
B2{1} = eval(B2j_R1);
C1{1} = eval(C1j_R1);

%                 A=B1^(-1)*A1
%                 B=B1^(-1)*B2
%                 C=B1^(-1)*C1
i=1;
            A=B1{i}\A1{i};
            B=B1{i}\B2{i};
            C=B1{i}\C1{i};  

%% The MOD Method under Fixed Regime 
      % regime 1
            [DET,FCC,OmegaK,~,IRF]=fmlre(A,B,C);  
            if isnan(DET(1,1)), [DET,~,~,~,Geig]=qzmlre(A,B,C); end
 
               if DET(1,1)<1 && DET(1,2)<=1,                         
                         disp('     (1)  Solution is DET'),disp(DET)    
                         disp([ ' , since ', ...
                             num2str(DET(1,1)) '<1 and ' num2str(DET(1,2)) '<=1'     ]);   disp(' ');         end                    
                % (2) Area for INDET                    
                    if DET(1,1)<1 && DET(1,2)>1,                         
                         disp('     (2)  Solution is  INDET'),disp(DET)   ; 
                          disp(' ');          end
                % (3) Area for NSS                               
                    if DET(1,1) >1 ,                         
                         disp('     (3)  Solution is   NSS'),disp(DET) ;   
                             disp(' ');        end
                         
                         