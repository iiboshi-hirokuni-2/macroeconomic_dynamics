%
%  Chapter 1
%  Basic Solow Model
%  Sec 1.5 page 14- 

load('K.mat')
 exp_k=k;

% initial setting 
rng('default')
nperiod = 100;
eps = 0.2*randn(nperiod,1); 
k = zeros(nperiod,1);
y = zeros(nperiod,1);

% parameters setting
delta = 0.1;
sigma = 0.2;
n = 0.02;
theta = 0.36;

% page 16
B = (1+theta*n-delta*(1-theta))/(1+n);
disp([ 'B = ',  num2str(B)]);

C = (delta+n)/(1+n);
disp([ 'C = ',  num2str(C)]);

% eq 1.4
for t = 1:nperiod
  k(t+1) = B*k(t) + C*eps(t); 
end  

% output
for t = 1:nperiod
  y(t+1) = theta*k(t) + eps(t); 
end  

% plot
figure('Name','Dynamics of Log Caipal')
  plot(k, 'b-','LineWidth',2);
  hold on
    plot(y, 'r-.','LineWidth',2);
  hold off
  title('Dynamics of Log Caipal and Output')
  xlabel('time');
  xlim([1, nperiod]);
  ylabel('Capital, Output')
  legend({'k', 'y'})
  set(gca, 'FontSize',12)


% Fig 1.3
figure('Name','Dynamics of Caipal')
  %plot(kss*(1+k), 'r-.','LineWidth',2);
   plot(kss*exp(k), 'r-.','LineWidth',2);
  hold on
    plot(exp_k, 'b-','LineWidth',2);
  hold off
  title('Dynamics of Log Linear and Exact model')
  xlabel('time');
  xlim([1, nperiod]);
  ylabel('Capital, Output')
  legend({'log k', 'K'})
  set(gca, 'FontSize',12)

