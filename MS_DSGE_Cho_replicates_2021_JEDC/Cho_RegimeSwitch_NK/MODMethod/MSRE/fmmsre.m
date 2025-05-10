function [DETC,FCC,OmegaK,GammaK,FK,OtherOutput,IRS]=fmmsre(P,A,B,C,R,Opt) 
%% MAJOR Updates of the previous code, written by Seonghoon Cho, March 25, 2019
% - Original Forward Method (FM) for LRE models was based on Cho and Moreno (2011,JEDC)
% - FM for MSRE models was based on Cho(2016,RED)
%   : Informally modified such that the forward solution coincides with the
%   MOD solution in the class of Determinacy-Admissible Models.
% - This is formally a "Modified" FM based on Cho(2019,WP),
%   implementing the MOD method for MSRE model.
%   *** Original FM still serves as a solution refinement.
%   *** Can be downloaded from http://web.yonsei.ac.kr/sc719/
%% This code
%   (1) checks the existence of the forward solution OmegaK,
%       computes it and the corresponding GammaK and FK
%   (2) provides conditions DETC for Classification of MSRE models
%       : determinacy/indeterminacy/no stable solution  
%% Acronyms 
%   LRE / MSRE : Linear Rational Expectations / Markov-Switching Rational Expectations
%   MOD Solution / Method : Minimum of Modulus (The most stable) solution / Method 
%   FM / FS : Forward Method / Forward Solution
%   DETC : Determinacy Conditions computed by FS.
%   FCC : Forward Convergence Condition
%   DET/INDET/NSS : : determinacy/indeterminacy/no stable solution
%   DA / DIA : Determinacy-Admissible / Determinacy-Inadmissible
%% THIS IS THE MODIFIED FORWARD METHOD (FM).  
%   - Original FM is FM under partial information, and is now option. 
%   - FM under Full information is also an option for FM
%   - This Modified FM (Default) is a sequential implementation of the following
%       1) Apply FM under partial information. 
%          Let the forward solution(FS) be FS_0. 
%          Computes DETC=[DETC1 DETC2 DETC3 DETC4] using FS_0. 
%          If DETC1*DETC2<1, the model is DA, FS=FS_0. Done.
%       2) If DETC1*DETC2>=1, Apply FM under full information.
%          Let the forward solution(FS) be FS_1. 
%          Computes DETC=[DETC1 DETC2 DETC3 DETC4] using FS_1. 
%          Define FS= argmin(DETC1 under partial information, DETC1 under full information).
%% Sequential Implementation of the MOD method 
%% : Decision making about Classification of MSRE models.
%   [1] Forward Method using this code.
%       (1) If DETC1*DETC2<1, the model is DA. FS=MOD. Done.
%       (2) If DETC1*DETC2>=1 and DETC1<1, model is indeterminate. FS=MOD is not confirmed.
%       (3) If DETC1*DETC2>=1 and DETC1>=1, FM cannot tell DET/INDET/NSS.
%   [2] Apply Groebner Basis(GB) Approach using gbmsre.m in the case of (3) from [1]
%       : Identify the MOD solution by solving all MSV solutions, completing the MOD method. 
%   NOTE 1: So far, without exception
%           1) every single model in the case of [1]-(1) is DA.
%           2) every single model in the remaining cases is DIA(Determinacy-Inadmissible).
%           Therefore, FM is sufficient in practice to identify Determinacy.
%   Note 2: GB approach is a stand-alone complete technique for implementing the MOD method. 
%           However, it is extremely time-consuming, not working most of
%           models with dimension greater than 2 or regimes greater than 2.
%% References :
%   Cho and Moreno(2011,JEDC), "Forward Method as a Solution Refinement in Rational Expectations Models"
%   Cho and McCallum(2015,JMacro),"Refining Linear Rational Expectations Models and Equilibria"
%   Cho(2016,RED) "Sufficient Conditions for Determinacy in a Class of Markov-Switching Rational 
%         Expectations Models" and its " A Technical Guide"'
%   Cho(2019,WP) "Determinacy and Classification of Markov-Switching Rational Expectations Models"
%   Foerster, et.al.(2016,QE) "Perturbation methods for Markov-switching dynamic stochastic general equilibrium models,"
%% Usage: Specify the input arguments and run it
%       P,A (required) 
%       B,C,R and other optional arguments (optional).
%   The code will produce 
%       DETC : Determinacy information+condition for determinacy-admissible model 
%       FCC : forward convergence information
%       OmegaK,GammaK,FK: Forward Solution
%       Other output : Functions of OmegaK, FK, and Impulse Response Functions
%% Code Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   The Class of Markov-Switching Rational Expectations model:
%   B1(s(t))x(t)=E[A1(s(t),s(t+1))x(t+1)|I(t)]+B2(s(t))*x(t-1)+C1(s(t))*z(t)   
%           x(t)=E[A(s(t),s(t+1))x(t+1)|I(t)]+B(s(t))*x(t-1)+C(s(t))*z(t)
%           z(t)=R(s(t))*z(t-1)+e(t)
%           where x(t) is an n by 1 vector of endogenous variables
%                 z(t) is an m by 1 vector of exogenous variables
%                 s(t) is a S - state Markov chain.
%                 A(s(t),s(t+1))=B1(s(t))^(-1)*A1(s(t),s(t+1))
%                 B(s(t))=B1(s(t))^(-1)*B2(s(t))
%                 C(s(t))=B1(s(t))^(-1)*C1(s(t))
%   Full set of RE solutions 
%       x(t)=Omega(s(t))*x(t-1)+Gamma(s(t))*z(t)+w(t), 
%       w(t)=E[F(s(t),s(t+1)w(t+1)|I(t)]
%           When w(t)=0,  x(t) is an MSV solution. For each Omega, there is
%                         a unique corresponding F.
%           When w(t)~=0, x(t) is a sunspot(non-fundamental) solution.
%   Definitions: 
%       (1) Omega_MOD = min DETC1(Omega) for all Omega where r() is the spectral radius 
%           F_MOD is the corresponding F  
%       (2) The model is Determinacy-Admissible iff Omega_MOD is real-valued and DETC1(Omega_MOD)*DETC2(F_MOD)<1
%           Therefore, if DETC1(OmegaK)*DETC2(FK), FS=MOD and the model is DA.
%       (3) OmegaK = forward solution and FK is the corresponding F.
%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       [1-1] Required
%           P   : S by S transition matrix where the ij-th element, 
%                 p_ij=Pr(s(t+1)=j|s(t)=i).
%           A   : S by S cell array such that A{i,j}=A(s(t)=i,s(t+1)=j), or
%                 S by 1 cell array such that A{i,1}=A(s(t)=i)
%                 ** size(A{i,j}) = n by n
%       [1-2] Optional
%           B   : S by 1 cell array such that B{i,1}=B(s(t)=i). 
%                   - Default: B=zeros(n,n);
%                   - If there is no x(t-1) in the model, set B=[]; 
%           C   : S by 1 cell array such that C{i,1}=C(s(t)=i).
%                   - size(C{i,1})=n by m.
%                   - Default: C{i,1}=eye(n) for all i, so m=n; 
%                   - If C is not present, set C=[];
%                   - C(s(t))=B0(s(t))^(-1)*C1(s(t) IS STRONGLY RECOMMENDED
%                     : Refer to Technical note on how to write a model.
%           R   : m by m VAR coefficient matrix of z(t) or
%                 S by 1 cell array such that R{i,1}=R(s(t)=1). (NEW)
%                   - Default: R=zeros(m,m);
%                   - If R is zeros, set R=[]; 
%           Opt : Other option-structure
%               Opt.maxK: maximum number of forward recursion: 
%                   (Default:maxK=1000)
%                   Change this only if the code display a warning sign. 
%               Opt.tolK: Precision of Forward Convergence of Omega_k: 
%                   The forward iteration stops if the difference between 
%                   the largest element of Omega_k at all states in 
%                   absolute value at k and k-1 is less than Opt.tolK.
%                   (Default is 0.000001)
%               Opt.IRsigma : m by 1 (or 1 by m ) vector of standard errors
%                   of the structural shocks, e(t) used for generating IR 
%                   functions. (Default : m by 1 vector of ones)
%               Opt.IRT: number of period for IR functions. (Default=20)
%               Opt.IRvs: 
%                   Default arrangement : responses of the first variable 
%                   to m shocks, ...,responses of the last variable to m shocks. 
%                   Set Opt.IRvs=1 if users want to arrange IR function 
%                   such that the responses of n variables to the first 
%                   shock, ..., responses of n variables to the last shock. 
%         (NEW) Opt.InfoAdj:  Whether to adjust information structure
%                   Opt.InfoAdj=0 : Original Cho-Moreno (2011) FM (Partial Information)
%                   Opt.InfoAdj=1 : FM (full information) 
%                   Opt.InfoAdj=-1 Modified FM (Default) : Automatic implementation:
%                   SEE technical Note
%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%       DETC : [DETC1 DETC2 DETC3 DETC4]
%               DETC1=r(BarPsi_kron(OmegaK,OmegaK)) : 
%               DETC2=r(Psi_kron(FK,FK)).
%               DETC3=DETC1*DETC2
%               DETC4=max(abs(eig(Psi_kron(OmegaK',FK)));
%       FCC : [FCC1 FCC2]   
%               FCC1 = termK : number of forward iteration for Omega_K. 
%                   If termK < maxK, the FCC holds for Omega_K
%               FCC2 = r(Psi_kron(R',F_K)) : n*m*S by n*m*S  
%                   If R_Psi_RtkFK < 1, the FCC holds for Gamma_K
%       OmegaK  : S by 1 cell array such that OmegaK{i,1}=Omega_k(s(t)=i) 
%                   at "k=termK". If termK < maxK, OmegaK is the forward
%                   solution.
%                   If termK = maxK, OmegaK is not the forward solution.
%       GammaK  : S by 1 cell array such that GammaK{i,1}=Gamma_k(s(t)=i) 
%                   at "k=termK". If FCC2 < 1, GammaK is the forward
%                   solution associated with OmegaK.
%                   If FCC2>=1, z(t) is explosive, thus GammaK is not
%                   the forward solution.
%       FK      : S by S cell array such that FK{i,j}=F_k(s(t)=i,s(t+1)=j) 
%                   at "k=termK".
%       OtherOutput : Structure of Other output information 
%               OtherOutput.BarPsi_OKkOK : the matrix BarPsi_OKkOK
%               OtherOutput.Psi_FKkFK : the matrix Psi_FKkFK
%               OtherOutput.Psi_RtkFK : the matrix Psi_RtkFK
%               OtherOutput.Psi_OKtkFK : the matrix Psi_OKtkFK
%       IRS : S by 1 cell array where i-th element of IRS is the IR of x(t)
%              when s(t)=i. It is T by n*m where each column is 
%               x(1)....x(n) to 1st shock,...  x(1)....x(n) to m-th shock 
%% Interpretation of the Results
%   (1) If FCC1<maxK, FS exists.
%       (1-1) If DETC1*DETC2<1. Then, FS=MOD and the model is DA. 
%               - If DETC1<1 and DETC2<1, the model is DET.  
%               - If DETC1<1 and DETC2.=1, the model is INDDET.
%               - If DETC1>=1, the model has NSS.
%       (1-2) If DETC1*DETC2>=1 and DETC1<1. the model is INDET.
%       (1-3) If DETC1*DETC2>=1 and DETC1>=1, Classification is not possible by FM.
%   (2) If FCC1="NaN", forward iteration reaches maxK, FS does not exist.
%       Classification is not possible by FM.
%       Apply the GB approach (gbmsre.m).

%% [1] Handling Input
%% [1-1] Identify S (number of Regimes), n (rows of A) and m (columns of C)
    S=size(P,1); 
    n=size(A{1,1},1);
    if size(A,2)==1, AA=cell(S,S); 
        for i=1:S, for j=1:S, AA{i,j}=A{i,1}; end, end, A=AA; 
    end
  
    if nargin==2, B=[]; C=[]; R=[]; end
    if nargin==3,       C=[]; R=[]; end
    if nargin==4,             R=[]; end
    
%% [1-2] Make B and R cell arrays if they are not cell arrays. 
    if isempty(B),      for j=1:S, B{j,1}=zeros(n,n);       end, end  
    if isempty(C),      for j=1:S, C{j,1}=eye(n,n);         end, end  
                                    m=size(C{1,1},2);
    if isempty(R),                  R=zeros(m);             end  
    if iscell(R), RR=R;  
    else, RR=cell(S,1); for i=1:S, RR{i,1}=R; end 
    end, R=RR;
                                    
%% [1-3] Option values  
    if nargin==6 && isfield(Opt,'maxK'), maxK=Opt.maxK; else, maxK=1000; end
    if nargin==6 && isfield(Opt,'tolK') 
                        tolK=Opt.tolK; else, tolK=0.000001; end
    if nargin==6 && isfield(Opt,'IRsigma') 
                        vecsig=Opt.IRsigma; else, vecsig=ones(1,m); end   
    if nargin==6 && isfield(Opt,'IRT'), T=Opt.IRT; else, T=20; end
    if nargin==6 && isfield(Opt,'IRvs'), vs=Opt.IRvs; else, vs=0; end
    if nargin==6 && isfield(Opt,'InfoAdj'), InfoAdj=Opt.InfoAdj; else, InfoAdj=-1; end

    
%% [2] Forward Method 
%% [2-1] Extension of the Original Cho-Moreno (2011) Forward Method (FM) to MSRE models
    if InfoAdj==0, [DETC,FCC,OmegaK,GammaK,FK,OtherOutput]=fm0(A,B);  end
%% [2-2] FM under full information 
    if InfoAdj==1, [DETC,FCC,OmegaK,GammaK,FK,OtherOutput]=fm1(A,B); end  
%% [2-3] Modified Forward Method (Default)
    % This Algoritm tries FM0 First.
    % If DETC(3)=r(BarPsi_kron(OmegaK,OmegaK))*r(Psi_kron(FK,FK))<1 under FM0, stop.
    % If not, apply FM1
    %    If DETC(3)<1 under FM1, stop. 
    %    If not, report results with the one that minimizes DETC(3) between FM0 and FM1
    if InfoAdj==-1, [DETC,FCC,OmegaK,GammaK,FK,OtherOutput]=fm(A,B); end

%% [3] Impulse-Responses Function under Regime-Switching
if nargout>6    
    D12=diag(vecsig); 
        %FT{j,t) and GT{j,t} are the coefficients of x(t) and z(t) 
        %in E_0[(x(t)].
        % Construct FT{j,1) and GT{j,1)
                for i=1:S
                    FT{i,1}=zeros(n,n); GT{i,1}=zeros(n,m); 
                    for j=1:S                        
                        FT{i,1}=FT{i,1}+P(i,j)*OmegaK{j};
                        GT{i,1}=GT{i,1}+P(i,j)*GammaK{j}*RR{j}; 
                 end, end
      
            for t=2:T-1
                for i=1:S
                    FT{i,t}=zeros(n,n); GT{i,t}=zeros(n,m);
                    for j=1:S                      
                        FT{i,t}=FT{i,t}+P(i,j)*FT{j,t-1}*OmegaK{j};
                        GT{i,t}=GT{i,t}+P(i,j)*((GT{j,t-1}+FT{j,t-1}...
                                *GammaK{j})*R{j}); 
                end, end
            end
        
        for i=1:S
                if vs==1 
                    tmp=GammaK{i}*D12; IRs=tmp(:)';
                        for t=1:T-1, temp=(FT{i,t}*GammaK{i}+GT{i,t})*D12; 
                            IRs=[IRs;temp(:)']; 
                        end
                    IRS{i,1}=IRs;
                else % Default
                    tmp=(GammaK{i}*D12)'; IRs=tmp(:)';
                        for t=1:T-1
                            temp=((FT{i,t}*GammaK{i}+GT{i,t})*D12)'; 
                            IRs=[IRs;temp(:)']; 
                        end
                    IRS{i,1}=IRs;    
                end
        end               

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nested Functions 
% ======================================================================
% Function 1 : FM0: Original Forward Method
% ======================================================================
    function [DETC,FCC,OmegaK,GammaK,FK,OtherOutput]=fm0(A,B)
    FS_Exist=1;       
    % Initial 目(s(t)=i) at k=1
        for i=1:S, OMEGAK{1,i}=B{i};  end
    % Fill out Xi_k(s(t)),目_k(s(t)), 
        k=0; FC_O=1;
        while FC_O>tolK
            k=k+1;   
            % (1) Xi{k,i}=eye(n)-E[A(s(t),s(t+1))*目_{k-1}(s(t+1))| s(t)=i]         
                for i=1:S 
                    EAOmegaki=zeros(n,n); 
                        for j=1:S                
                           EAOmegaki =EAOmegaki +P(i,j)*A{i,j}*OMEGAK{k,j};
                        end
                    Xi{k,i}=eye(n)-EAOmegaki;
                end

     % Compute 目_{k}(s(t)) 
                FC_O_OmegaS=[]; 
                for i=1:S 
                    OMEGAK{k+1,i}=Xi{k,i}\B{i};
                    FC_O_OmegaS=[FC_O_OmegaS;...
                                max(abs(OMEGAK{k+1,i}(:)-OMEGAK{k,i}(:)))];
                end
            FC_O=max([FC_O_OmegaS]);     
            % Checking whether FS is finite
            if k==maxK, FS_Exist=0; break; end         
            for s=1:S
                if max(abs(OMEGAK{k,s}(:)))>10e10, FS_Exist=0; break; end
            end
        end % end of "while".
        
        termK=k;   
        OmegaK=cell(S,1); for i=1:S, OmegaK{i,1}=OMEGAK{end,i};  end  
    % Computing All of the Results from OmegaK at k=termK
        [DETC,FCC,OmegaK,GammaK,FK,OtherOutput]=fm_results(OmegaK,termK,FS_Exist);
    end

% ======================================================================
% Function 2 : FM1: Full Information Forward Method
% ======================================================================
    %      Make agents use all variables of the the model as state variables 
    %      when forming expectations. 
    function [DETC,FCC,OmegaK,GammaK,FK,OtherOutput]=fm1(A,B)     
        tt=100; % Just in case the FCC holds, but with a large k such that termK=maxK
              % this step may be done more than once.
        for t1=1:tt
            if t1<tt/2, H=0.1*t1*ones(n,n); % 
            else, H=10*randn(n,n); %              
            end
            [OmegaK,termK,FS_Exist]=fm1_OmegaK(A,B,H); 
            [DETC]=fm_results(OmegaK,termK,FS_Exist);
            if FS_Exist==1 && DETC(3)<1, break, end  % OmegaK is the MOD.     
        end       
        [DETC,FCC,OmegaK,GammaK,FK,OtherOutput]=fm_results(OmegaK,termK,FS_Exist);
    end

    %    Computing OmegaK under Full Information: to be read by fm1.m.
    
    function [OmegaK,termK,FS_Exist]=fm1_OmegaK(A,B,H)
        FS_Exist=1; 
        
        EB=cell(S,1);
        for i=1:S
            EBi=zeros(n,n); 
            for j=1:S              
                EBi =EBi +P(i,j)*B{j,1};
            end
            EB{i,1}=EBi;
        end
        
        for i=1:S
            BigB1{i,1}=[eye(n) -eye(n);H*EB{i,1} eye(n)];
            BigB2{i,1}=[B{i,1} zeros(n,n);zeros(n,2*n)];
            for j=1:S          
                BigA1{i,j}=[zeros(n,2*n);A{i,j}+H -H];
                BigA{i,j}=BigB1{i,1}\BigA1{i,j}; 
            end
            BigB{i,1}=BigB1{i,1}\BigB2{i,1}; 
        end    
        N=2*n;

    % Initial 目(s(t)=i) at k=1
        for i=1:S, BigOMEGAK{1,i}=BigB{i,1};  end
    % Fill out Xi_k(s(t)),目_k(s(t)), 
        k=0; FC_O=1;
        while FC_O>tolK
            k=k+1;   
            % (1) Xi{k,i}=eye(n)-E[A(s(t),s(t+1))*目_{k-1}(s(t+1))| s(t)=i]         
                for i=1:S
                    EAOmegaki=zeros(N,N); 
                        for j=1:S               
                           EAOmegaki =EAOmegaki +P(i,j)*BigA{i,j}*BigOMEGAK{k,j};
                        end
                    BigXi{k,i}=eye(N)-EAOmegaki;
                end

            % (4) Compute 目_{k}(s(t)) 
                FC_O_OmegaS=[]; 
                for i=1:S 
                    BigOMEGAK{k+1,i}=BigXi{k,i}\BigB{i};
                    FC_O_OmegaS=[FC_O_OmegaS;...
                                max(abs(BigOMEGAK{k+1,i}(:)-BigOMEGAK{k,i}(:)))];
                end
            FC_O=max([FC_O_OmegaS]);     
            
            % Checking whether FS is finite
            if k==maxK, FS_Exist=0; break; end         
            for s=1:S
                if max(abs(BigOMEGAK{k,s}(:)))>10e10, FS_Exist=0; break; end
            end 
            
        end % end of "while".
        
        termK=k;   
        OmegaK=cell(S,1); for i=1:S, OmegaK{i,1}=BigOMEGAK{end,i}(1:n,1:n);  end   
        
    end

% ======================================================================
% Function 3 : Modified Forward Method
% ======================================================================
    function [DETC,FCC,OmegaK,GammaK,FK,OtherOutput]=fm(A,B) % 
    % Let FM0 (FM1) be Original (Full info) FM.
    % (1) Try FM0 first.
        [DETC_0,FCC_0,OmegaK_0,GammaK_0,FK_0,OtherOutput_0]=fm0(A,B); 
        % (1-1) If DETC_0(3)<1 under FM0, Done.
        if DETC_0(3)<1, DETC=DETC_0; FCC=FCC_0; OmegaK=OmegaK_0; GammaK=GammaK_0; FK=FK_0; OtherOutput=OtherOutput_0;
        % (1-2-1) If NaN, Apply FM1, Done        
        elseif isnan(DETC_0(3)), [DETC,FCC,OmegaK,GammaK,FK,OtherOutput]=fm1(A,B); 
        % (1-2-2) If DETC(3)>=1 Apply FM1      
        elseif DETC_0(3)>=1, [DETC_1,FCC_1,OmegaK_1,GammaK_1,FK_1,OtherOutput_1]=fm1(A,B);
            % (1-2-2-1) Choose the Results that has smaller DETC(1)
            %           under FM0 and FM1.   
            if DETC_1(3)<1;              DETC=DETC_1; FCC=FCC_1; OmegaK=OmegaK_1; 
                                         GammaK=GammaK_1; FK=FK_1;OtherOutput=OtherOutput_1;
            elseif isnan(DETC_1(3)),     DETC=DETC_0; FCC=FCC_0; OmegaK=OmegaK_0; 
                                         GammaK=GammaK_0; FK=FK_0;OtherOutput=OtherOutput_0;
            elseif DETC_1(3)>=1
                if DETC_0(1)<=DETC_1(1), DETC=DETC_0; FCC=FCC_0; OmegaK=OmegaK_0; 
                                         GammaK=GammaK_0; FK=FK_0;OtherOutput=OtherOutput_0;
                else,                    DETC=DETC_1; FCC=FCC_1; OmegaK=OmegaK_1; 
                                         GammaK=GammaK_1; FK=FK_1;OtherOutput=OtherOutput_1;
                end
            end
        end
    end

% ======================================================================
% Function 4 : Common Outputs from the computed OmegaK
% ======================================================================
    function [DETC,FCC,OmegaK,GammaK,FK,OtherOutput]=fm_results(OmegaK,termK,FS_Exist)      
        if FS_Exist==1
             for i=1:S
                EAOmegaK=zeros(n,n); 
                    for j=1:S               
                        EAOmegaK =EAOmegaK +P(i,j)*A{i,j}*OmegaK{j,1};
                    end
                    XiK{i,1}=eye(n)-EAOmegaK;
            end 
            for i=1:S, for j=1:S, FK{i,j}=XiK{i,1}\A{i,j}; end, end

            % The remaining results
            BarPsi_OKkOK=BarPsi_X1X1(OmegaK);
            Psi_FKkFK=Psi_XX(FK);
            Psi_RtkFK=Psi_X1tY(R,FK);
            Psi_OKtkFK=Psi_X1tY(OmegaK,FK);
        
            DETC1=max(abs(eig(BarPsi_OKkOK)));
            DETC2=max(abs(eig(Psi_FKkFK)));
            MODWD=DETC1*DETC2;
            DETC4=max(abs(eig(Psi_OKtkFK)));
            FCC1=termK;  
            FCC2=max(abs(eig(Psi_RtkFK)));
               
            DETC=[DETC1 DETC2 MODWD DETC4];
            FCC=[FCC1 FCC2];
            
            % Gamma
            vvC=[]; 
            for i=1:S
                InvXiC=XiK{i,1}\C{i,1};
                vvC=[vvC;InvXiC(:)]; 
            end
            vvGamma=(eye(S*n*m)-Psi_RtkFK)\vvC;
            vGamma=reshape(vvGamma,n*m,S); GammaK=cell(S,1);
            for i=1:S, GammaK{i,1}=reshape(vGamma(:,i),n,m); end
            
        elseif FS_Exist==0
            DETC=[NaN NaN NaN NaN]; FCC=[NaN NaN]; 
            for i=1:S
                OmegaK{i,1}=NaN*ones(n,n); 
                GammaK{i,1}=NaN*ones(n,m); 
                for j=1:S, FK{i,j}=NaN*ones(n,n); end
                BarPsi_OKkOK=NaN*ones(n^2*S,n^2*S);
                Psi_FKkFK=NaN*ones(n^2*S,n^2*S);
                Psi_RtkFK=NaN*ones(n*m*S,n*m*S);
                Psi_OKtkFK=NaN*ones(n^2*S,n^2*S);
            end
        end
        OtherOutput.BarPsi_OKkOK=BarPsi_OKkOK;
        OtherOutput.Psi_FKkFK=Psi_FKkFK;
        OtherOutput.Psi_RtkFK=Psi_RtkFK;
        OtherOutput.Psi_OKtkFK=Psi_OKtkFK;
    end

% ======================================================================
% Other Functions : Kronecker Products of Matrices
% ======================================================================
    %% BarPsi_kron(OmegaK,OmegaK)
    function BarPsi_XkX=BarPsi_X1X1(X)
        bdiagXX=zeros(n^2*S,n^2*S);
        for i=1:S
            bdiagXX(n^2*(i-1)+1:n^2*i,n^2*(i-1)+1:n^2*i)...
            = kron(X{i},X{i});
        end
        BarPsi_XkX= bdiagXX*kron(P',eye(n^2));  %
    end

    %% Psi_kron(FK,FK) 
    function Psi_XkX=Psi_XX(X)
        Psi_XkX=[]; 
        for i=1:S
            Psi_XkXrow=[]; 
            for j=1:S
                Psi_XkXrow=[Psi_XkXrow P(i,j)*kron(X{i,j},X{i,j})];
            end
            Psi_XkX=[Psi_XkX;Psi_XkXrow]; %Psi_(kron(X,X)) 
        end
    end
    %% Psi_kron(R',FK) 
    function Psi_X1tkY=Psi_X1tY(X1,Y)
        Psi_X1tkY=[]; 
        for i=1:S
            Psi_X1tkYrow=[]; 
            for j=1:S
                    Psi_X1tkYrow=[Psi_X1tkYrow P(i,j)*kron(X1{j,1}',Y{i,j})];
            end
            Psi_X1tkY=[Psi_X1tkY;Psi_X1tkYrow];  
        end
    end
% ======================================================================
% end of Nested Function
% ======================================================================
end  %% end of main function