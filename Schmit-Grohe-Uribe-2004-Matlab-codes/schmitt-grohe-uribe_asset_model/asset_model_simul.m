             
%% simulation
% cal   unconditional_mean   
sig = 1;   % St Dev of shock
[Ey,Ex] = unconditional_mean(gx, hx, gxx, hxx, gss, hss, eta, sig)

% cal stochastic steady state by fixed point 
x_d = zeros(1,1);  % x_d = zeros(4*2,1);
 x0_ss = csolve('eq_stochastic_steady_state',x_d ,[],1e-4,1000,hx, hxx, hss, sig )
%   x0_ss2   = x0_ss(1:4,1)+x0_ss(5:8,1)
[Error] = eq_stochastic_steady_state(x0_ss , hx, hxx, hss, sig )

% Impulse Response Function
sig = 100;          % St Dev of shock
x0= 0;        % initial values of predetermined variables
h = 8;        % horizon of IRF
e = [ 1; zeros(h,1)]; % generations of shocks

[Y2,X2] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0_ss, e);
[Y1,X1] = simu_1st(gx, hx, eta, sig, x0_ss, e);
t = 0:1:h+1;

figure('Name', 'IRF of Asset Model')
subplot(2,1,1)
     plot(t,Y2-Y2(1,1),'r-','LineWidth',1.5);
   hold on
     plot(t,Y1,'b:','LineWidth',1.5);
   hold off  
     title('Control Variable: Price-Dividend Ratio, p_t/d_t')
     legend('2nd order', '1st order')
     xlim([0 h])
subplot(2,1,2)    
      plot(t,X2(:,1),'r-','LineWidth',1.5);
   hold on
     plot(t,X1(:,1),'b:','LineWidth',1.5);
   hold off  
     title('Pre-determined Variable: Growth Rate of Dividends, x_t')
     legend('2nd order', '1st order')
      xlim([0 h])


% simulation
sig = 100;          % St Dev of shock
x0= [Ex];        % initial values of predetermined variables
e = randn(100,1); % generations of shocks
% e = [ 1; zeros(19,1)];

[Y2,X2] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e);
[Y1,X1] = simu_1st(gx, hx, eta, sig, x0, e);

figure('Name', 'Artificial Data of Asset Model')
subplot(2,1,1)
     plot(Y2-Y2(1,1),'r-','LineWidth',1.5);
   hold on
     plot(Y1,'b:','LineWidth',1.5);
   hold off  
     title('Control Variable: Price-Dividend Ratio, p_t/d_t')
     legend('2nd order', '1st order')
     xlim([1 100])
subplot(2,1,2)    
      plot(X2(:,1),'r-','LineWidth',1.5);
   hold on
     plot(X1(:,1),'b:','LineWidth',1.5);
   hold off  
     title('Pre-determined Variable: Growth Rate of Dividends, x_t')
     legend('2nd order', '1st order')
      xlim([1 100])
      
   
      
      