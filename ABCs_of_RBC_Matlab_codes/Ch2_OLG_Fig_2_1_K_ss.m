%
% Chapter 2
% OLG
% Sec 2.2

clear all


% parameters setting
kappa = 4;
theta = 0.36;

% Page 27:  steady state of capital
i =1; j=1;

for k = 0:0.05:10
  k1 = kappa*k^theta;  
  K(i,:) = [k, k1];
  i = i+1;
  if (k1 < k)&(j==1)
     kss = k; 
     disp([ 'k_ss = ', num2str(k)]);
     j=0;
  end    
end

% Plot Figure 2.1
figure('Name','Fig 2.1')
  plot(K(:,1), K(:,2),'b-','LineWidth',2);
  hold on
    plot(K(:,1), K(:,1),'r-.','LineWidth',2);    
    xline(kss,'LineWidth',1.25);
    yline(kss,'LineWidth',1.25);
  hold off  
  title('Fig 2.1');
  xlabel('k_t');
  ylabel('k_{t+1}')
  legend({'Capital','45 degree', 'Steady State'},...
       'Location','northwest');
  set(gca, 'FontSize',12)
