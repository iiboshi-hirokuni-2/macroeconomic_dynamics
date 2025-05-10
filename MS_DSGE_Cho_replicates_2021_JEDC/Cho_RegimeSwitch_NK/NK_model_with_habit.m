% endogenous variables
syms  y ygap ynat rstar pai mc R b tau cons em g a nu dy_obs dc_obs pi_obs rn_obs
% future variables
syms  yf ygapf ynatf rstarf paif mcf Rf bf tauf consf emf  gf af nuf dy_obsf dc_obsf pi_obsf rn_obsf
% lag variables
syms  y1 ygap1 ynat1 rstar1 pai1 mc1 R1 b1 tau1 cons1  em1  g1 a1 nu1 dy_obs1 dc_obs1 pi_obs1 rn_obs1
% exogenous variables
syms eps_mp  eps_g eps_a eps_nu

% parameters
syms  kap kappa beta beta1 sig hc varphi rhoR psipi psiy b_st g_st phi_y
syms  rhog rhoa rhotau rhomp rhonu
syms  deltab deltay0  deltaestar R_st
%% make MATRIX
% Y_f =  transpose([ yf paif Rf bf  tauf  ef     ]);
% Y_t = transpose([ y pai R b tau e   ]);  
% Y_b = transpose([ y1 pai1 R1 b1  tau1  e1   ]);

Y_f =  transpose([ yf ygapf ynatf rstarf paif mcf Rf bf  tauf consf  emf  gf af nuf dy_obsf dc_obsf pi_obsf rn_obsf]);
Y_t = transpose([ y ygap ynat rstar pai mc R b tau cons em g a nu dy_obs dc_obs pi_obs rn_obs]);  
Y_b = transpose([ y1 ygap1 ynat1 rstar1 pai1 mc1 R1 b1 tau1 cons1 em1  g1 a1 nu1 dy_obs1 dc_obs1 pi_obs1 rn_obs1]);
Epsilon_t = transpose([  eps_mp eps_g eps_a eps_nu]);

% Model
% [1] Linearized Euler Equation
eq1 = - (1+hc)*ygap + ygapf + hc*ygap1- (1/sig)*(1-hc)*( R - paif - rstar);
%eq1 = - ygap + ygapf  - (1/sig)*( R - paif - rstar);

% [2] Phillips Curve  
eq2 = -pai + kap*mc + beta*paif + nu;

% [3] Marginal Cost 
eq3 = -mc + (kappa + (sig*hc/(1-hc)))*ygap -  (sig*hc/(1-hc))*ygap1;

% [4] Monetary policy rule 
eq4 = -R + (psiy)*(1-rhoR)*y+psipi*(1-rhoR)*pai + rhoR*R1 + em;

% [5] Gov Debt --> Eq5
% eq6 = -b + betta1*b1+ b_st*betta1*( - y + y1  -pai)- tau + e  ;
eq5 = -b + beta1*b1+ b_st*beta1*(R1 - y + y1  -pai)- tau + 1/g_st*g;

% [6] Tax Rule  --> Eq7
% eq5 = -tau + rhotau*tau1 +deltab*(1-rhotau)*b1 + deltay0*(1-rhotau)*y+ eps_tau;
eq6 = -tau + rhotau*tau1 +deltab*(1-rhotau)*b1 + deltay0*(1-rhotau)*ygap;

% [7] Output gap
eq7 = - ygap + y - ynat;

% [8] G shock
 eq8 = -g + rhog*g1 +eps_g; 

% [9] Technology shock
 eq9 = -a +   rhoa*a1 + eps_a;

% [10] NKPC shock
 eq10 = -nu +   rhonu*nu1 + eps_nu; 

% [11] MP shock
 eq11 = -em +   rhomp*em1 + eps_mp; 

% [12] Natural output
 eq12 = -(kappa + (sig*hc/(1-hc)))*ynat + (sig*hc/(1-hc))*ynat1 + (1+1/varphi)*a + (sig/(1-hc))*(g-hc*g1);

% [13] Natural rate 
 eq13 = -(1-hc)*rstar + sig*((1-rhog+hc)*g - hc*g1) - sig*((1+hc)*ynat - ynatf - hc*ynat1);

% [14] Cosumption
eq14 = -cons + y - g;

% [15-18] Obs Equations
eq15 = -dy_obs + y - y1;
eq16 = -dc_obs + cons - cons1;
eq17 = -pi_obs + pai;
eq18 = -rn_obs + R;
 


% Set Sims Forms
System_of_Eq = [  ...
                eq1;  eq2;  eq3;  eq4;  eq5;...
                eq6;  eq7;  eq8;  eq9;  eq10;  ...
                eq11; eq12; eq13; eq14; ...
                eq15; eq16; eq17; eq18; ...
                ];


neq  = length(System_of_Eq);
CC = zeros(neq,1);
% chek if all f evaluated at ss value are zero
check = eval(System_of_Eq)