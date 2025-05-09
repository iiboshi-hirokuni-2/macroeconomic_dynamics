
%% Check Solution   


% cal   unconditional_mean   
sig = 3;          % St Dev of shock
[Ey,Ex] = unconditional_mean(gx, hx, gxx, hxx, gss, hss, eta, sig)


% stochastic steady state by iteration
x0 = Ex;
e = [ 0 ; zeros(200,1)];
[Y02,X02] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e);
x_sss = X02(end,:)' % stochastic steady state 

% cal stochastic steady state by fixed point 
 x0_ss = csolve('eq_stochastic_steady_state',Ex,[],1e-10,50,hx, hxx, hss, sig )
 [Error] = eq_stochastic_steady_state(x0_ss , hx, hxx, hss, sig )

% Impulse Response Function
% sig = 1;     % St Dev of shock
x0_1= [0;0];      % initial values of predetermined variables
% x0_2 =Ex; 
x0_2 = x0_ss;
h = 10;        % horizon of IRF
e = [ 1; zeros(h,1)]; % generations of TFP shocks
t = 0:1:h+1;

[Y2,X2] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0_2, e);
[Y1,X1] = simu_1st(gx, hx, eta, sig, x0_1, e);

figure('Name','IRF of Neoclassical Model')
subplot(3,1,1)
     plot(t,Y2-Y2(1,1),'r-','LineWidth',1.5);
   hold on
     plot(t,Y1-Y1(1,1),'b:','LineWidth',1.5);
   hold off  
     title('Control Variable: log c_t')
     legend('2nd order', '1st order')
     xlim([0 h])
%      ylim([-0.1 2])
subplot(3,1,2)    
      plot(t,X2(:,1)-X2(1,1),'r-','LineWidth',1.5);
   hold on
     plot(t,X1(:,1)-X1(1,1),'b:','LineWidth',1.5);
   hold off  
     title('Pre-determined Variable: log k_t')
     legend('2nd order', '1st order')
      xlim([0 h])
%       ylim([-0.1 2])
subplot(3,1,3)    
      plot(t,X2(:,2),'r-','LineWidth',1.5);
   hold on
     plot(t,X1(:,2),'b:','LineWidth',1.5);
   hold off  
     title('Pre-determined Variable: log A_t')
     legend('2nd order', '1st order')
      xlim([0 h])

% Generation of Artificial Data
% sig = 1;          % St Dev of shock
x0 =[kp; ap];    % initial values of predetermined variables
                 % kp and ap are steady states
period = 100;                 
e = randn(period,1); % generations of shocks


[Y2,X2] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e);
[Y1,X1] = simu_1st(gx, hx, eta, sig, x0, e);

figure('Name','Artificial Data of Neoclassical Model')
subplot(3,1,1)
     plot(Y2,'r-','LineWidth',1.5);
   hold on
     plot(Y1,'b:','LineWidth',1.5);
   hold off  
     title('Control Variable: log c_t')
     legend('2nd order', '1st order')
     xlim([1 period])
subplot(3,1,2)    
      plot(X2(:,1),'r-','LineWidth',1.5);
   hold on
     plot(X1(:,1),'b:','LineWidth',1.5);
   hold off  
     title('Pre-determined Variable: log k_t')
     legend('2nd order', '1st order')
      xlim([1 period])
subplot(3,1,3)    
      plot(X2(:,2),'r-','LineWidth',1.5);
   hold on
     plot(X1(:,2),'b:','LineWidth',1.5);
   hold off  
     title('Pre-determined Variable: log A_t')
     legend('2nd order', '1st order')
      xlim([1 period])
     
     
% print anal_deriv
filename =('neoclassical_model');
ETASHOCK = 1;

anal_deriv_print2f(filename,fx,fxp,fy,fyp,f,ETASHOCK, ...
    fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx)
   


