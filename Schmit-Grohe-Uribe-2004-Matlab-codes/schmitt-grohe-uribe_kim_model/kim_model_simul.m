%% Check Solution    
% cal   unconditional_mean   
sig = 1;  % St Dev of shock
[Ey,Ex] = unconditional_mean(gx, hx, gxx, hxx, gss, hss, eta, sig)

x0 = Ex;
e = [ 0 0; zeros(200,2)];
[Y02,X02] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e);
x_ss1 = X02(end,:)' % stochastic steady state 

% cal stochastic steady state by fixed point 
x_d = zeros(4,1);  % x_d = zeros(4*2,1);
 x0_ss = csolve('eq_stochastic_steady_state',x_d ,[],1e-5,2000,hx, hxx, hss, sig )
%   x0_ss2   = x0_ss(1:4,1)+x0_ss(5:8,1)
[Error] = eq_stochastic_steady_state(x0_ss , hx, hxx, hss, sig )


% Impulse Response Function
% sig = 0.1;          % St Dev of shock
x0_1 = [0;0;0;0];        % initial values of predetermined variables
x0_2 = x_ss1; 
h = 40;        % horizon of IRF
e = [ 1 0; zeros(h,2)]; % generations of shocks
t = 0:1:h+1; 

[Y2,X2] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0_2, e);
[Y1,X1] = simu_1st(gx, hx, eta, sig, x0_1, e);

figure('Name','IRF of Kim Model')
subplot(3,1,1)
     plot(t,Y2-Y2(1,1),'r-','LineWidth',1.5);
   hold on
     plot(t,Y1-Y1(1,1),'b:','LineWidth',1.5);
   hold off  
     title('Control Variable: log c_t')
     legend('2nd order', '1st order')
     xlim([0 h])
%      ylim( [-0.025 0.025])
subplot(3,1,2)    
      plot(t,X2(:,1)-X2(1,1),'r-','LineWidth',1.5);
   hold on
     plot(t,X1(:,1)-X1(1,1),'b:','LineWidth',1.5);
   hold off  
     title('Pre-determined Variable: log k_t')
     legend('2nd order', '1st order')
      xlim([0 h])
%       ylim( [-1 2])
subplot(3,1,3)    
      plot(t,X2(:,3),'r-','LineWidth',1.5);
   hold on
     plot(t,X1(:,3),'b:','LineWidth',1.5);
   hold off  
     title('Pre-determined Variable: log A_t')
     legend('2nd order', '1st order')
      xlim([0 h])


% simulation
% sig = 0.5;          % St Dev of shock
x0= zeros(4,1);        % initial values of predetermined variables
period = 100;
e = randn(period,2); % generations of shocks

[Y2,X2] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e);
[Y1,X1] = simu_1st(gx, hx, eta, sig, x0, e);

figure('Name','Artificial Data of Kim model')
subplot(3,1,1)
     plot(Y2-Y2(1,1),'r-','LineWidth',1.5);
   hold on
     plot(Y1-Y1(1,1),'b:','LineWidth',1.5);
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
      plot(X2(:,3),'r-','LineWidth',1.5);
   hold on
     plot(X1(:,3),'b:','LineWidth',1.5);
   hold off  
     title('Pre-determined Variable: log A_t')
     legend('2nd order', '1st order')
      xlim([1 period])

