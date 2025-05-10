%% This example illustrates computation time of Groebner basis (GB) technique
%% Summary of Results from Experiments
%   - Singular is based on C++ and known as more than 100 times faster than matlab. 
%   - Nevertheless, GB produces solves the problem (identifying the MOD
%     solution by computing all MSV solutions) within a very limited set of models.
%   - For model(n,m,S) with a dimension n, m lagged endogenous variables and S regimes.
%     GB works within an hour in the following cases 
%            LRE Models(n,m,1)    MSRE Models(n,m,2)   MSRE Models(n,m,3)
%               n<24, m=1            n<=6, m=1          n<=3,  m=1 
%               n<=6, m=2            n=2,  m=2                      
%               n<=4, m=3            
%   - When increasing any one of m, n and S, computation stops, 
%     out of memory after around a day using a standard PC.
%   - This experiments shows the reason.
%           Def: N=n*m*S be the number of unknowns of the set of 
%                the second order polynomials implied by a model(n,m,S)
%           Def: nsol denotes the number of solutions. 
%       1) Increase in N by 1 increases computational time by a factor of 2^nsol.
%       2) A rise in n, m or S increases nsol by change in n*m*S. 
%       3) Thus, nsol increases computational time by a factor of 2^(n*m*S)
%   - Therefore, GB would be very difficult to apply for standard macroeconomic models.
%% Idea of the experiments 
%   - Imagine that in LRE models, we have no information about the eigensystem, so QZ method is unknown.
%   - This is the restriction we face when solving for MSRE models.
%   - This enables us to gauze the approximate time increase by increasing N. 
%% Writing the Problem in Polynomials
% 	(1) LRE model 
%       The problem can be written as a set of second order polynomials.
%       Let x(t)=Omega*x(t-1) be an MSV solution to a model x(t)=AE[x(t+1)]+Bx(t-1).
%       The problem is to solve for Omega such that A*Omega^2-Omega+B=0
%       Vectorizing this yields a set of n^2 polynomials in n^2 unknowns.
%       Let m=rank(B). Then there are m*n unknowns in total, including
%       m second order unknowns(non-zero diagonal elements of Omega)
%   (2) MSRE model      
%       The problem can be written as a set of second order polynomials.
%       Let x(t)=Omega(s(t)*x(t-1) be an MSV solution to an MSRE model x(t)=E[A(s(t),s(t+1))x(t+1)]+B(s(t)x(t-1).
%       The problem is to solve for Omega(i) 
%           (Sum_j P(i,j)A(i,j)Omega(j))*Omega(i)-Omega(i)+B(i)=0 for all i=1,..S.
%       Vectorizing this yields a set of S*n^2 polynomials in S*n^2 unknowns.
%       Let m=rank(B(s(t)). Then there are m*n*S unknowns in total, including
%       m*S second order unknowns(non-zero diagonal elements of Omega(s(t)))
%% Number of Solutions and Difficulties
%   (1) LRE model(n,m) 
%       - The qz method shows that maximum nsol= C(2n-(n-m),n)=C(n+m,n),
%         choosing n out of n+m.
%   (2) MSRE model(n,m,S)
%       - nsol is finite but unknown. This can be interpreted
%         as a special case of S*n^2 unkowns in as many polynomials with
%         some restrictions.
%       - Therefore, maximum nsol < C(2nS,nS), but propotional to this.
%       - For instance, 
%                              nsol , C(2nS,nS)   Computation time
%           n=1, m=2, S=2,       4        6        < a second
%           n=2, m=2, S=2,      44       70        about 227 seconds
%           n=3, m=2, S=2,      ??      924        Never stops     
%% Instruction
%% Choose model by setting, n,m and S. Set S=1 for LRE, 2 or higher for MSRE
%   For LRE models,  set n=1,2,3,.., m=1,2,3,...,S=1 and run this code.
%   For MSRE models, set n=1,2,3,.., m=1,2,3,.., S=2,3, and run this code.
%   * You may also type [DETC,FCC,OmegaK]=fmmsre(P,A,B); to see the result from FM.

    clear
    n=2; m=1; S=2;
%   This code generates a random MSRE model of the specified dimensions. 
%   The LRE model is the one with the first regime being permanent. 
%   See below.

%% Summary of the Results for LRE models 
%% Computation time increases by a factor of 2^nsol as nsol increases
%    Model  n   m    S    #(sols)   max nsol(m=n)  GB works?  #(ideal) Time to compute
%     LRE   1   1    1          2           2        Yes           1     Within a second
%     LRE   2   1    1          3           6        Yes           3     Within a second
%     LRE   3   1    1          4          20        Yes           6     Within a second
%     LRE   4   1    1          5          70        Yes          10     Within a second
%     LRE   5   1    1          6         252        Yes          15     Within a second  
%     LRE   6   1    1          7         924        Yes          21     Within a second  
%     LRE   7   1    1          8        3432        Yes          28     Within a second  
%     LRE   8   1    1          9       12800        Yes          36     Within a second  
%     LRE   9   1    1          10      48620        Yes         45     Around 1.5 second  
%     LRE   10   1    1         11     184756        Yes         55     Around 3 second 
%     LRE   11   1    1         12     705432        Yes         66     Around 4 second 
%     LRE   12   1    1         13    2704156        Yes         78     12 seconds 
%     LRE   13   1    1         14   10400600        Yes         91     35 seconds
%     LRE   14   1    1         15   40116600        Yes         105    95 seconds
%     LRE   15   1    1         16  155117520        Yes         120    208 seconds
%     ..................................................................................
%     LER   24   1    1         25  3.2248e+13        No     24*25/2    More than a day expected.  
%     ==========================================================================================
%     LRE   2   2    1          6           6        Yes          10     Within a second
%     LRE   3   2    1          10         20        Yes          21     Within a second
%     LRE   4   2    1          15         70        Yes          37     1 second 
%     LRE   5   2    1          21        252        Yes          59     7 seconds
%     LRE   6   2    1          28        924        Yes          86     548 seconds, increase by a factor of 2^(n+1)    
%     LRE   7   2    1          36       3432        No           ??     38 hours expected increase by a factor of 2^(n+1)      
%     ==========================================================================================
%     LRE   3   3    1          20         20        Yes          47     3 seconds. 
%     LRE   4   3    1          35         70        Yes          87     More than an hour
%     LRE   5   3    1          56        252        No          ??      Not within a day
%     ==========================================================================================
%     LRE   4   4    1          70         70        No           ??     Not within a day

%% Summary of the Results for MSRE models 
%% Computation time increases by a factor of 2^nsol as nsol increases
%    Model  n   m    S    #(sols)   max nsol(m=n)  GB works?  #(ideal) Time to compute                             
%     MSRE  1   1    2          4        <=6        Yes            3     Within a second
%     MSRE  2   1    2          9        <=70       Yes           11     Within a second
%     MSRE  3   1    2         16       <=924       Yes           26     Within a second
%     MSRE  4   1    2         25     <=12800       Yes           50     3 seconds
%     MSRE  5   1    2         36    <=184756       Yes           85     60 seconds
%     MSRE  6   1    2         49   <=2704156       Yes           133    573 seconds
%     MSRE  7   1    2         64  <=40116600        No           196?   ?? More than a day (2048 minutes) expected
%     ==========================================================================================
%     MSRE  2   2    2          44       <=70        Yes          61     227 seconds
%     MSRE  3   2    2         >=100     <=924       No           ??     Never within a day
%     ==========================================================================================
%     MSRE  3   3    2         >=400      < ??       No           ??     Never within a day
%     ==========================================================================================
%     MSRE  1   1    3          8        <=20        Yes           7     Within a second
%     MSRE  2   1    3          27       <=924       Yes          32     2 seconds
%     MSRE  3   1    3          64     <=48620       Yes          93     592 seconds
%     MSRE  4   1    3          ??    <=48620       Yes           ??     ??
%     ==========================================================================================
%     MSRE  2   2    3         ??       <=924        No           ??     Not within a day

%% Random model
    if S==1 
        disp(strjoin(["LRE Model,", "n=", num2str(n), "m=", num2str(m)]))
    elseif S>1
        disp(strjoin(["MSRE Model,", "n=", num2str(n), "m=", num2str(m), "S=", num2str(S)]))
    end


    P=0.1*rand(S,S); for Si=1:S, P(Si,Si)=0.9+P(Si,Si); end
    P=P./(P*ones(S,S));
    for Si=1:S, for Sj=1:S, A{Si,Sj}=0.5*randn(n,n); end, end
    for Si=1:S, B{Si,1}=0.4*randn(n,n); end
    

    
    if m<n, for Si=1:S, B{Si,1}(:,1:n-m)=zeros(n,n-m); end, end

    [DETCMOD,OmegaMOD,FMOD,DETC_All,AllOmegas,DETC_OmjtkF1]=gbmsre(P,A,B);  
