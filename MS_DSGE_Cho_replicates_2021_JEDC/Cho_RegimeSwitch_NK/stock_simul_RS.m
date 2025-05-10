function [State_path,Endo_path,Exo_path] = stock_simul_RS(OmegaK,GammaK,P,init_state,T, vecsig)

n_stata = max(size(OmegaK));
n_endo = max(size(OmegaK{1}));
n_exo = min(size(GammaK{1}));
D12=diag(vecsig);

sum_P = tril(ones(size(P)))*P;
State_path = zeros(1,T);
Endo_path = zeros(n_endo,T);
Exo_path = zeros(n_exo,T);


State_path(1) = init_state;
for t = 2:T
    state0 = State_path(t-1);    
    n = rand();
    state1 = find(sum_P(:,state0)-n>0,1);
    State_path(t) = state1;
    
    Shock = diag(D12.*randn(n_exo));
    Exo_path(:,t)= Shock;
    
    Endo_past = Endo_path(:,t-1);    
    Endo_current = OmegaK{state1}*Endo_past + GammaK{state1}*Shock;
    Endo_path(:,t) = Endo_current;
end


