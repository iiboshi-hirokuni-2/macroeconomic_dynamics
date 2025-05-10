 a0=0.98; % Prob(S_t=1|S_t-1=1)
 b1=0.975; % Prob(S_t=2|S_t-1=2)
 
 PP = [
        a0      1-b1
        1-a0     b1
     ];