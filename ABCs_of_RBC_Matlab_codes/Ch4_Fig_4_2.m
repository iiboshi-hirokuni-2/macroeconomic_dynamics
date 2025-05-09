%
%  Chapter 4
%
%


global vlast beta delta theta k0 kt
hold off
hold all

%set initial conditions
vlast=zeros(1,100);
k0=0.06:0.06:6;
beta=.98;
delta=.1;
theta=.36;
numits=240;

%begin the recursive calculations
for k=1:numits
    for j=1:100
        kt=j*.06;
        %find the maximum of the value function (minus the minimum)
        ktp1=fminbnd(@valfun,0.01,6.2);
        v(j)=-valfun(ktp1);
        kt1(j)=ktp1; % policy function
    end
    if k/48==round(k/48)
        %plot the steps in finding the value function
        figure(100)
        hold on
          plot(k0,v)
          yline(0)
        hold off
        drawnow
        title('Fig 4.2 value function')
        xlabel('k_t');
        ylabel('V(k_{t})')
        set(gca,'FontSize',14)
    end
    vlast=v;
end


%% plot the final policy function
for i = 1:length(k0)
   if kt1(i) < k0(i)
      kss = k0(i);
      break
   end
end

figure('Name','Fig 4.3 policy function')
plot(k0,kt1,'LineWidth',2)
hold on
  plot(k0,k0,'LineWidth',2)  
  xline(kss,'LineWidth',1.5) ; 
hold off
title('Fig 4.3 policy function')
xlabel('k_t');
ylabel('k_{t+1}')
set(gca,'FontSize',14)

%% function
function val=valfun(k)

global vlast beta delta theta k0 kt
%smooth out the previous value function
g=interp1(k0,vlast,k,'linear');
%Calculate consumption with given paraameters
kk=kt^theta-k+(1-delta)*kt;
if kk <= 0
    %to keep values from going negative
    val=-888-800*abs(kk);
else
    %calculate the value of the value function at k
    val=log(kk)+beta*g;
end
%change value to negative since program finds minimum
val=-val;

end



