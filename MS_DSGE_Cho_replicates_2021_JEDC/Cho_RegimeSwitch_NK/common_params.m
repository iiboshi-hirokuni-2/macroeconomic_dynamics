beta        =0.99;
beta1       = 1/beta;
sig         = 1.5;
%sig         = 1.0;
theta       = 0.75;
varphi      = 1/2.5;
%varphi      = 1/1;
hc          = 0.7;   %habit
%alpha_0     = 1/3;
%ave_mat     = 20;

kappa       = 1/varphi + sig;
%Cap_phi     = 0.8616;
kap         = (1-theta)*(1-theta*beta)/theta;
%deltay0     = 0.2823;

rhog        = 0.95;   
rhoa        = 0.50;
rhomp       = 0.1;
rhonu       = 0.75;
    

pi_st      = 0;
gam_st      = 0;
%b_st        = 1.25;
b_st        = 1.25;
g_st        = 1/(1-0.20);
%tau_st      = 0.171379820612903;    

% Steady state real interest rate
%rr_steady   = (1+gam_st) / beta-1;
rr_steady   = beta1-1;
% Steady state FFR (gross)
%R_st        = 1+(pi_st + rr_steady);
R_st        = 1 + rr_steady;

    
%iota_y      = 0;
%phi_y       = 0.99;

%rhowood     = (1-1/ave_mat)*betta1;

%if (rhowood>beta1)
%    Error; 
%end

% Stds of structural shocks
p_er = [ 0.0019;
    0.0028;
    0.0062;
    0.0030];

% Stds of observation errors
o_er = zeros(3,1);




