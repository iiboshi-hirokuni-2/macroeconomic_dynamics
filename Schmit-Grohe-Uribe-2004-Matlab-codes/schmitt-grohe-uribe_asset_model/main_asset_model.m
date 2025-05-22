close all
clear all

%ASSET_RUN.M
%[gx,hx,gxx,hxx,gss,hss] = asset_run_linear(eta) computes a first- or second-order approximation (depending on the value assigned to the variable approx in the program ASSET_MODEL.M)to the policy function  of the Burnside asset pricing model as described in ``Accuracy of stochastic perturbation methods: the case of asset pricing models," by F. Collard and M. Juillard, JEDC, 25, June 2001, 979-999. 
%
%(c) Stephanie Schmitt-Grohe and Martin Uribe
%
%January 22, 2002
%
%
%The solution is of the form
%
%xp = h(x,sigma) + sigma* eta * ep
%
%y = g(x,sigma)
%
%The quadratic approximation to these functions are (in tensor notation) [Notation: x is x_t and xp is x_t+1, variables are expressed in log-deviations from their steady state value]
%
%xp^i = hx^i_a x_a + 1/2 [hxx^i_ab x_a x_b + hss^i sigma^2] + sigma*  eta^i_c ep_c
%
%y^i = gx^i_a x_a + 1/2 [gxx^i_ab x_a x_b + gss^i sigma^2]
%
%where 
%
% hx is nx by nx
% gx is ny by nx
% hxx is nx by nx by nx
% gxx is ny by nx by nx
% eta is nx by ne
% gss is ny by 1
% hss is nx by 1
% sigma is a positive scalar
%
%Calls: asset_model.m gx_hx.m gxx_hxx.m gss_hss.m
% 
%(c) Stephanie Schmitt-Grohe and Martin Uribe
%Date January 22, 2002

addpath('./function/')

P.RHO=-0.139; P.THETA=-1.5; P.SIG =1; 

%  analytical 2nd order derivatives of system: f
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,...
    fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,...
    fxy,fxxp,fxx,eta]=asset_model(P); 

%First-order approximation of Policy function
[gx,hx] = gx_hx(fy,fx,fyp,fxp)

%Second-order approximation of Policy function
[gxx,hxx] = gxx_hxx(fx,fxp,fy,fyp,fypyp,fypy,fypxp,...
                   fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,...
                   fxpxp,fxpx,fxyp,fxy,fxxp,fxx,hx,gx) 

[gss,hss] = gss_hss(fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,...
                    fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,...
                    fxy,fxxp,fxx,hx,gx,gxx,eta)

%% Simulation of solution of model
asset_model_simul;

%% Check Table 1 for 2nd-order approximation 
disp(' ')
P.THETA=-1.5; P.RHO=-0.139;  P.SIG =1; 
asset_model_table_1;
%
disp(' ')
P.THETA=-10; P.RHO=-0.139; 
asset_model_table_1;
                
%
disp(' ')
P.THETA=-1.5; P.RHO=0.9; 
asset_model_table_1;   

%% stochastic steady state by iteration
x0 = Ex;
e = [ 0 ; zeros(200,1)];
[Y02,X02] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e);
x_ss = X02(end,:)'; % stochastic steady state 

%% cal stochastic steady state by fixed point 
x_d = zeros(1,1);  % x_d = zeros(4*2,1);
 x0_ss = csolve('eq_stochastic_steady_state',x_d ,[],1e-4,1000,hx, hxx, hss, sig );
%   x0_ss2   = x0_ss(1:4,1)+x0_ss(5:8,1)
[Error] = eq_stochastic_steady_state(x0_ss , hx, hxx, hss, sig );
