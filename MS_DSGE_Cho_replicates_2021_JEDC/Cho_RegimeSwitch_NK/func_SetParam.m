
function [A,B,C]=func_SetParam(par,B1j,A1j,B2j,C1j )

%=====================
common_params;

% regime 1 (PMP/AFP) 
Regime_parameters_R1;

%========================
psipi      =  par(1);
deltab     = par(3);
%========================

B1{1} = eval(B1j);
A1{1} = eval(A1j);
B2{1} = eval(B2j);
C1{1} = eval(C1j);

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
psipi       =  par(2);
deltab     = par(4);
%========================

B1{2} = eval(B1j);
A1{2} = eval(A1j);
B2{2} = eval(B2j);
C1{2} = eval(C1j);

%                 A=B1^(-1)*A1
%                 B=B1^(-1)*B2
%                 C=B1^(-1)*C1
i=2;
            A{i,1}=B1{i}\A1{i};
            B{i,1}=B1{i}\B2{i};
            C{i,1}=B1{i}\C1{i};  
 








