% endogenous variables
syms  y ygap ynat rstar pai R b tau em g a
% future variables
syms  yf ygapf ynatf rstarf paif Rf bf  tauf  emf  gf af
% lag variables
syms y1 ygap1 ynat1 rstar1 pai1 R1 b1 tau1  em1  g1 a1
% exogenous variables
syms eps_mp  eps_g eps_a

% parameters
syms  kap kappa beta beta1 sig varphi rhoR psipi psiy b_st g_st phi_y
syms  rhog rhoa rhotau rhomp 
syms  deltab deltay0  deltaestar R_st
%% make MATRIX
% Y_f =  transpose([ yf paif Rf bf  tauf  ef     ]);
% Y_t = transpose([ y pai R b tau e   ]);  
% Y_b = transpose([ y1 pai1 R1 b1  tau1  e1   ]);

Y_f =  transpose([ yf ygapf ynatf rstarf paif Rf bf  tauf  emf  gf af ]);
Y_t = transpose([ y ygap ynat rstar pai R b tau em g a ]);  
Y_b = transpose([ y1 ygap1 ynat1 rstar1 pai1 R1 b1 tau1  em1  g1 a1 ]);
Epsilon_t = transpose([  eps_mp eps_g eps_a ]);

% Model
% [1] Linearized Euler Equation
eq1 = - ygap + ygapf  - (1/sig)*( R - paif - rstar);

% [2] Phillips Curve  
eq2 = -pai + kap*kappa*ygap + beta*paif;

% [3] Monetary policy rule 
eq3 = -R + (psiy)*(1-rhoR)*y+psipi*(1-rhoR)*pai + rhoR*R1 + em;

% [4] Gov Debt --> Eq5
% eq6 = -b + betta1*b1+ b_st*betta1*( - y + y1  -pai)- tau + e  ;
eq4 = -b + beta1*b1+ b_st*beta1*(R1 - y + y1  -pai)- tau + 1/g_st*g;

% [5] Tax Rule  --> Eq7
% eq5 = -tau + rhotau*tau1 +deltab*(1-rhotau)*b1 + deltay0*(1-rhotau)*y+ eps_tau;
eq5 = -tau + rhotau*tau1 +deltab*(1-rhotau)*b1 + deltay0*(1-rhotau)*ygap;

% [6] Output gap
eq6 = - ygap + y - ynat;

% [7] G shock
 eq7 = -g + rhog*g1 +eps_g; 

% [8] Technology shock
 eq8 = -a +   rhoa*a1 + eps_a;

% [9] MP shock
 eq9 = -em +   rhomp*em1 + eps_mp; 

% [10] Natural output
 eq10 = -ynat +((1+1/varphi)/kappa)*a + (sig/kappa)*g;

% [11] Natural rate 
 eq11 = -rstar + sig*(1-rhog)*g  - sig*(ynat - ynatf);


% Set Sims Forms
System_of_Eq = [  ...
                eq1;  eq2;  eq3;  eq4;  eq5;...
                eq6;  eq7;  eq8;  eq9;  eq10;  eq11; ...
                ];


neq  = length(System_of_Eq);
CC = zeros(neq,1);
% chek if all f evaluated at ss value are zero
check = eval(System_of_Eq)