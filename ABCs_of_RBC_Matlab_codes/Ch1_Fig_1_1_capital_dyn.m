%
%  Chapter 1
%  Basic Solow Model
%  Sec 1.1 - Sec 1.4 
%  page 9 - 12 

% initial setting 
rng('default')
nperiod = 100;
eps = 0.2*randn(nperiod,1); 
ngrid = 100;

% parameters setting
delta = 0.1;
sigma = 0.2;
n = 0.02;
theta = 0.36;
A_bar =1;

% eq 1.1: steady state of capital
i =1; j=1;

for k = 0:0.01:3
  k1 = ((1-delta)*k + sigma*A_bar*k^theta) /(1+n);  
  K(i,:) = [k, k1];
  i = i+1;
  if (k1 < k)&(j==1)
     kss = k; 
     disp([ 'k_ss = ', num2str(k)]);
     j=0;
  end    
end

% Plot Figure 1.1
figure('Name','Fig 1.1')
  plot(K(:,1), K(:,2),'b-','LineWidth',2);
  hold on
    plot(K(:,1), K(:,1),'r-.','LineWidth',2);    
    xline(kss,'LineWidth',1.25);
    yline(kss,'LineWidth',1.25);
  hold off  
  title('Fig 1.1');
  xlabel('k_t');
  ylabel('k_{t+1}')
  legend({'Capital','45 degree', 'Steady State'},...
       'Location','northwest');
  set(gca, 'FontSize',12)


% Sec 1.4 
% kss = 2.23;
k = ones(nperiod,1)*kss;
y = zeros(nperiod,1);

% eq 1.2
for t = 1:nperiod
  fk= k(t)^theta;
  k(t+1)= ((1-delta)*k(t)+sigma*A_bar*exp(eps(t))*fk)/(1+n);

  y(t+1) = A_bar*exp(eps(t))*fk;
end  

% plot
figure('Name','Dynamics of Caipal')
  plot(k, 'b-','LineWidth',2);
  hold on
    plot(y, 'r-.','LineWidth',2);
  hold off
  title('Dynamics of Caipal and Output')
  xlabel('time');
  xlim([1, nperiod]);
  ylabel('Capital, Output')
  legend({'k', 'y'})
  set(gca, 'FontSize',12)

  save('K.mat','k','kss')
