
nperiod = 200;
state = zeros(nperiod); 


p1=.8;
p2=1-p1;

for t = 1:nperiod
    if p1 > rand(1)
       state(t) =1 ;
    else
       state(t) =2 ;
    end
end    

figure(100)
plot(state)
ylim([0.5 2.5])
title('State of Tech')
set(gca, 'FontSize',12)

%% Markov chain

p11 = 0.8;
p21 = 1-p11;
p22 = 0.9;
p12 = 1-p22;

nperiod = 200;

state_mc = zeros(nperiod);
state_mc(1) = 1;

for t = 2:nperiod
    if state_mc(t-1)== 1
       if p11 > rand(1)
           state_mc(t) =1 ;
       else
           state_mc(t) =2 ;
       end
    else % state = 2
       if p21 > rand(1)
           state_mc(t) =1 ;
       else
           state_mc(t) =2 ;
       end 
    end   
end    

figure(200)
plot(state_mc)
ylim([0.5 2.5])
title('State of Tech by MC')
set(gca, 'FontSize',12)

