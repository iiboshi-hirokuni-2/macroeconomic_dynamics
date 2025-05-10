%% This code compares the computation time for
%% 1) qzmlre, 2) qzmlre2gensys.m (gensys) 3) fmlre for LRE models and 4) fmmsre for MSRE models.
%% 1 through N-dimensional MSRE model with S regime is randomly generated.
%% The first three codes take regime 1 as an LRE model.
%% Computation time is the average of T trials.
%% Results
%% 1) Whenever forward iteration is less than 50, forward method is
%%    faster than standard methods (QZ, gensys). 
%% 2) Computation time for a MSRE model several times longer than LRE counterpart.
%% 3) Computation time mildly increases as the dimension of the model is larger.


S=2;
N=10;
T=10;
TtableT=zeros(N,6);

for t=1:T
    Ttable=[];
for n=1:N
    P=0.1*rand(S,S); for Si=1:S, P(Si,Si)=0.9+P(Si,Si); end
    P=P./(P*ones(S,S));
    
    k=0; tmax=1; 
    while tmax>0.5
        k=k+1;
        for i=1:S
            for j=1:S
                A{i,j}=diag(0.5+0.5*rand(n,1))+0.05*randn(n,n); 
            end
            B{i,1}=diag(0.9*ones(n,1))-A{i,1}+0.05*randn(n,n); 
        end
                tic, DETCMOD=qzmlre(A{1,1},B{1,1}); t1=toc;  
                tic, [G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=qzmlre2gensys(A{1,1},B{1,1});  t2=toc;
                tic, [DETC1,FCC1]=fmlre(A{1,1},B{1,1});  t3=toc;
                tic, [DETC,FCC]=fmmsre(P,A,B);  t4=toc;
                
                if isnan(FCC1(1)), t3K=1; else, t3K=0; end
                if isnan(FCC(1)), t4K=1; else, t4K=0; end
                tmax=max([t3 t4 t3K t4K]);
    end
    
    Ttable=[Ttable;[t1 t2 t3 t4 FCC1(1) FCC(1)]]; 
end

TtableT=TtableT+Ttable;
t
end

TtableAVG=[(1:1:N)' TtableT/T];
disp('Computation time Comparison')
disp('   Dimension  qzmlre   gensys     fmlre     fmmsre     K_LRE    K_MSRE') 
disp(TtableAVG) 
disp('K_LRE and K_MSRE are numbers of forward iteration in Forward Method for LRE and MSRE models.')
    