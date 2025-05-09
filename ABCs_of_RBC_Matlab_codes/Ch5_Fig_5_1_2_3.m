%
%  Chapter 4
%
%

clear all
close all


global vlast1 vlast2 beta delta theta k0 kt At p1 p2
hold off
hold all

% initial value function
vlast1=20*ones(1,40);
vlast2=vlast1;

% initial policy function
k0=0.4:0.4:16;
kt11=k0;
kt12=k0;

% setting parameters
beta=.98;
delta=.1;
theta=.36;
A1=1.75;
A2=.75;
p1=.8;
p2=1-p1;


% calculate value function
numits=250;
for k=1:numits
    for j=1:40
        kt=k0(j);

      % Case A1
        At=A1;
        z=fminbnd(@valfunsto,.41,15.99);
        v1(j)=-valfunsto(z);
        kt11(j)=z; % policy function

      % Case A2
        At=A2;
        z=fminbnd(@valfunsto,.41,15.99);
        v2(j)=-valfunsto(z);
        kt12(j)=z; % policy function
    end

    if k/50==round(k/50)        
        figure(500)
        hold on
           plot(k0,v1,'b-')
           plot(k0,v2,'r-.')
          yline(0)
        hold off
        drawnow
        title('Fig 5.1 value function')
        xlabel('k_t');
        ylabel('V(k_{t})')
        ylim([25 50])
        set(gca,'FontSize',14)       
        drawnow
    end
    vlast1=v1;
    vlast2=v2;
end
hold off

%% value function
figure(500)
        hold on
           l1=plot(k0,v1,'b-','LineWidth',2);
           l2=plot(k0,v2,'r-.','LineWidth',2);
          yline(0)
        hold off
        drawnow
        title('Fig 5.1: value function')
        xlabel('k_t');
        ylabel('V(k_{t})')
        legend([l1,l2],{'A1','A2'},'Location','northwest');
        ylim([25 50])
        set(gca,'FontSize',14)

%% policy function

figure('Name','fig 5.2')
plot(k0,kt11,k0,kt12,'LineWidth',2)
title('Fig 5.2: policy function')
legend({'A1','A2'},'Location','northwest');
xlabel('k_t')
ylabel('k_{t+1}')
set(gca,'FontSize',14)

%% State of Technology: A
nperiod = 500;
state = zeros(nperiod); 


for t = 1:nperiod
    if p1 > rand(1)
       state(t) =1 ;
    else
       state(t) =2 ;
    end
end    

figure('Name','Trace of State')
plot(state)
ylim([0.75 2.25])
title('State of Technology: A')
xlabel('time')
ylabel('A_t')
set(gca, 'FontSize',12)

%% Fig 5.3: time path of capital

k_dyn = zeros(nperiod,1);
k_dyn(1) = 10;

for t = 2:nperiod
    k = k_dyn(t-1); 
    
    if state(t)==1 
      k_dyn(t) =  interp1(k0,kt11,k,'linear');
    else
       k_dyn(t) =  interp1(k0,kt12,k,'linear');
    end     
end   

figure('Name','Fig 5.3')
plot(k_dyn,'LineWidth',2)
% ylim([0.5 2.5])
xlabel('time')
ylabel('k_t')
title('Fig 5.3: time path of capital')
set(gca, 'FontSize',12)


%% function
function val=valfunsto(x)
global vlast1 vlast2 beta delta theta k0 kt At p1 p2
k=x;
g1=interp1(k0,vlast1,k,'linear');
g2=interp1(k0,vlast2,k,'linear');
kk=At*kt^theta-k+(1-delta)*kt;
if kk<=.001
    val=log(.001)+beta*(p1*g1+p2*g2)+200*(kk-.001);
else
    val=log(kk)+beta*(p1*g1+p2*g2);
end
val=-val;

end
