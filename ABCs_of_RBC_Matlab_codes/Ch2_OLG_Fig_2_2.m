%
% Chapter 2
% OLG
%

clear all 

% parameters setting
kappa = 4;
theta = 0.36;
roe = 0.9;

% initial capital
K0 = [0.8*kappa^(1/(1-theta)),...
      kappa^(1/(1-theta)), ...
      1.2*kappa^(1/(1-theta))];

lambda(1:3,1)=[1 1 1]';    

H =65;
K(1:3,1)=K0;
Y(:,1)=lambda(:,1).*K(:,1).^theta.*H.^(1-theta);


%  dynamics of OLG
for t=2:120
    shock = 0.02.*(rand(3,1)-0.5);
    lambda(:,t)=(1-roe)+roe.*lambda(:,t-1)+shock;

    % capital
    K(:,t)=kappa.*lambda(:,t).*K(:,t-1).^theta;

    % output
    Y(:,t)=lambda(:,t).*K(:,t-1).^theta.*H.^(1-theta);
end    

% Plot
figure('Name','Ch2')
subplot(2,1,1)
  plot(K','LineWidth',2)
  title('K')
  legend({ '0.8','1.0','1.2' });
  set(gca,'FontSize',12)
subplot(2,1,2)
  plot(Y','LineWidth',2)
  title('Y')
  legend({ '0.8','1.0','1.2' });
  set(gca,'FontSize',12)
