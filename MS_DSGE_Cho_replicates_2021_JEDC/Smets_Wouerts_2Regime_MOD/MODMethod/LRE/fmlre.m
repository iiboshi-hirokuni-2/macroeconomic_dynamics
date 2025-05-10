function [DETC,FCC,OmegaK,GammaK,FK,IRs]=fmlre(A,B,C,R,Opt)
%% MAJOR Updates of the previous code, written by Seonghoon Cho, March 25, 2019
% - Original Forward Method (FM) for LRE models was based on Cho and Moreno (2011,JEDC)
% - Informally modified FM is proposed following Cho(2016,RED) such that 
%   the forward solution coincides with the MOD solution in the class of Determinacy-Admissible Models.
% - This is formally a "Modified" FM based on Cho(2019,WP),
%   implementing the MOD method for LRE model.
%   *** Original FM still serves as a solution refinement.
%   *** Can be downloaded from http://web.yonsei.ac.kr/sc719/
%% This code
%   (1) checks the existence of the forward solution OmegaK,
%       computes it and the corresponding GammaK and FK
%   (2) provides conditions DETC for Classification of MSRE models
%       : determinacy/indeterminacy/no stable solution  
%% Acronyms 
%   LRE : Linear Rational Expectations 
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
%          Computes DETC=[DETC1 DETC2 DETC3] using FS_0. 
%          If DETC1*DETC2<1, the model is DA, FS=FS_0. Done.
%       2) If DETC1*DETC2>=1, Apply FM under full information.
%          Let the forward solution(FS) be FS_1. 
%          Computes DETC=[DETC1 DETC2 DETC3] using FS_1. 
%          Define FS= argmin(DETC1 under partial information, DETC1 under full information).
%% Sequential Implementation of the MOD method 
%% : Decision making about Classification of MSRE models.
%   [1] Forward Method using this code.
%       (1) If DETC1*DETC2<1, the model is DA. FS=MOD. Done.
%       (2) If DETC1*DETC2>=1 and DETC1<1, model is indeterminate. FS=MOD is not confirmed.
%       (3) If DETC1*DETC2>=1 and DETC1>=1, FM cannot tell DET/INDET/NSS.
%       (4) If DETC1=NaN, the forward fail to converge. FM cannot tell DET/INDET/NSS.
%   [2] Apply QZ method using qzmlre.m in the case of (3) from [1]
%       : Identify the MOD solution by solving all MSV solutions, completing the MOD method. 
%   NOTE 1: So far, without an exception
%           1) every single model in the case of [1]-(1) is DA.
%           2) every single model in the remaining cases is DIA(Determinacy-Inadmissible).
%           Therefore, FM is sufficient in practice to identify Determinacy.
%   Note 2: QZ method is a stand-alone complete technique for implementing the MOD method. 
%           QZ method may be preferred to this hybrid approach if the
%           research interests are the classification of the model. But FM
%           can be used as a equilibrium refinement as proposed by Cho and
%           Moreno (2011, JEDC)
%% References :
%   Cho and Moreno(2011,JEDC), "Forward Method as a Solution Refinement in Rational Expectations Models"
%   Cho and McCallum(2015,JMacro),"Refininglinearrationalexpectationsmodelsandequilibria"
%   Cho(2016,RED) "Sufficient Conditions for Determinacy in a Class of Markov-Switching Rational 
%         Expectations Models" and its " A Technical Guide"'
%   Cho(2019,WP) "Determinacy and Classification of Markov-Switching Rational Expectations Models"
%% Usage: Specify the input arguments and run it
%       A (required) 
%       B,C,R and other optional arguments (optional).
%   The code will produce 
%       DETC : Determinacy information+condition for determinacy-admissible model 
%       FCC : forward convergence information
%       OmegaK,GammaK,FK: Forward Solution
%       Other output : Impulse Response Functions
%% Code Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   The Class of Linear Rational Expectations model:
%         B1x(t)=A1*E[x(t+1)|I(t)]+B2*x(t-1)+C1*z(t)
%           x(t)=A*E[x(t+1)|I(t)]+B*x(t-1)+C*z(t)
%           z(t)=R*z(t-1)+e(t)
%           where x(t) is an n by 1 vector of endogenous variables
%                 z(t) is an m by 1 vector of exogenous variables
%                 A=B1^(-1)*A1
%                 B=B1^(-1)*B2
%                 C=B1^(-1)*C1
%   Full set of the RE solutions:
%       x(t)=Omega*x(t-1)+Gamma*z(t) +w(t), 
%           w(t)=F* E_t [w(t+1)]
%           When w(t)=0,  x(t) is an MSV solution. For each Omega, there is
%                         a unique corresponding F.
%           When w(t)~=0, x(t) is a sunspot(non-fundamental) solution.
%   Definitions: 
%       (1) Omega_MOD = min r(Omega) for all Omega where r() is the spectral radius 
%           F_MOD is the corresponding F  
%       (2) The model is Determinacy-Admissible iff Omega_MOD is real-valued and r(Omega_MOD)*r(F_MOD)<1
%       (3) OmegaK = forward solution and FK is the corresponding F.
%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       [1-1] Required 
%           A : n by n matrix (forward looking, may well be singular)
%       [1-2] Optional 
%           B : n by n matrix (backwardlooking, may well be singular)
%                           Default : B=zeros(n,n);
%           C : n by m coefficient matrix of the m by 1 vector of
%                exogenous variables, z(t). (Default: C=eye(n))
%           R :  VAR coefficient matrix of z(t-1) (Default:R=zeros(m,m))
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
%                   Default arrangement : first variable to m shocks, ..., 
%                   the last variable to m shocks. 
%                   Set Opt.IRvs=1 if users want to arrange IR function 
%                   such that the n variables to the first shock, ..., n 
%                   variables to the last shock. 
%         (New) Opt.InfoAdj:  Whether to adjust information structure
%                   Opt.InfoAdj=0 : Original Cho-Moreno (2011) FM (Partial Information)
%                   Opt.InfoAdj=1 : FM (full information) 
%                   Opt.InfoAdj=-1 Modified FM (Default) : Automatic implementation:
%                   SEE technical Note
%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         DETC :     [DETC1 DETC2 DETC1*DETC2]
%                   = [r(OmegaK) r(FK) r(OmegaK)*r(FK)] where r() is spectral radius.  
%         FCC :     [FCC1 FCC2] = [termK R_RtKFK]  
%                       termK : number of forward iteration for Omega_K. 
%                       If FCC1 < maxK, the FCC holds for Omega_K
%                       R_RtKFK= r(kron(R',FK)) : 
%                       If FCC2 < 1, the FCC holds for Gamma_k
%         OmegaK:   Omega_k at "k=termK". 
%                   If termK = maxK, OmegaK is not the forward solution.
%         GammaK:   Gamma_k at "k=termK". If FCC2 < 1, GammaK is Gamma^*
%                   If FCC2>=1, GammaK is not the forward solution.
%         FK:       F_k at "k=termK".
%         IRS :     IR function of x(t)
%                   It is T by n*m where each column is 
%                   x(1)....x(n) to 1st shock,...  x(1)....x(n) to m-th shock 
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
%       Apply the QZ method (qzmlre.m).

%% [1] Handling Input
%% [1-1] Identify n (rows of A) and m (columns of C)
    n=size(A,1);
    if nargin==1, B=[]; C=[]; R=[]; end
    if nargin==2,       C=[]; R=[]; end
    if nargin==3,             R=[]; end
%% [1-2] Define B, C or R if they are not given. 
    if isempty(B),      B=zeros(n,n); end  
    if isempty(C),      C=eye(n,n);   end  
                        m=size(C,2);
    if isempty(R),      R=zeros(m,m); end  
    

%% [1-3] Option values  
    if nargin==5 && isfield(Opt,'maxK'), maxK=Opt.maxK; else maxK=1000; end
    if nargin==5 && isfield(Opt,'tolK') 
                        tolK=Opt.tolK; else, tolK=0.000001; end
    if nargin==5 && isfield(Opt,'IRsigma') 
                        vecsig=Opt.IRsigma; else, vecsig=ones(1,m); end   
    if nargin==5 && isfield(Opt,'IRT'), T=Opt.IRT; else, T=20; end
    if nargin==5 && isfield(Opt,'IRvs'), vs=Opt.IRvs; else, vs=0; end
    if nargin==5 && isfield(Opt,'InfoAdj'), InfoAdj=Opt.InfoAdj; else, InfoAdj=-1; end
   
    
%% [2] Forward Method

%% [2-1] FM0: Original Cho (2016) Forward Method (FM) for MSRE Models
    if InfoAdj==0, [DETC,FCC,OmegaK,GammaK,FK]=fm0(A,B); end
    
%% [2-2] FM1: FM under full information 
    if InfoAdj==1, [DETC,FCC,OmegaK,GammaK,FK]=fm1(A,B); end 
    
%% [2-3] Modified Forward Method (Default)
    % This Algoritm tries FM0 First.
    % If DETC(3)=r(OmegaK)*r(FK)<1 under FM0, stop.
    % If not, apply FM1
    %    If DETC(3)<1 under FM1, stop. 
    %    If not, report results with the one that minimizes DETC(3) between FM0 and FM1
    if InfoAdj==-1, [DETC,FCC,OmegaK,GammaK,FK]=fm(A,B); end

%% [3] Other Output: Impulse-Responses Functions 
    if nargout>5    
        D12=diag(vecsig); 
        % FT(t) and GT(t) are the coefficients of x(t) and z(t) in E_0[(x(t)].
          FT{1,1}=zeros(n,n); FT{1,1}=OmegaK; 
          GT{1,1}=zeros(n,m); GT{1,1}=GammaK*R;                                     
            for t=2:T-1
               	FT{1,t}=zeros(n,n); FT{1,t}=FT{1,t-1}*OmegaK;     
                GT{1,t}=zeros(n,m); GT{1,t}=(GT{1,t-1}+FT{1,t-1}*GammaK)*R; 
            end
            
            if vs==1
               tmp=GammaK*D12; IRs=tmp(:)';
                for t=1:T-1
                    temp=(FT{1,t}*GammaK+GT{1,t})*D12; 
                    IRs=[IRs;temp(:)']; 
                end
            else % Default
                tmp=(GammaK*D12)'; IRs=tmp(:)';
                for t=1:T-1 
                    temp=((FT{1,t}*GammaK+GT{1,t})*D12)'; 
                    IRs=[IRs;temp(:)']; 
                end
            end           
    end

%% End of the Main Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nested Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% ======================================================================
% Function 1 : FM0: Original Forward Method
% ======================================================================
    function [DETC,FCC,OmegaK,GammaK,FK]=fm0(A,B) % 
        FS_Exist=1;    
        % Initial 目_k at k=1
            OMEGAK{1,1}=B;
        % 目_k for all k 
            k=0; FC_O=1; 
            while FC_O>tolK
                k=k+1;   
                    OMEGAK{k+1,1}=(eye(n)-A*OMEGAK{k,1})\B;
                    FC_O=max(abs(OMEGAK{k+1,1}(:)-OMEGAK{k,1}(:)));     
                % Checking whether FS is finite    
                if k==maxK, FS_Exist=0; break; end     
                if max(abs(OMEGAK{k,1}(:)))>10e10, FS_Exist=0; break; end
            end 
        termK=k; OmegaK=OMEGAK{end,1};
        
        [DETC,FCC,OmegaK,GammaK,FK]=fm_results(OmegaK,termK,FS_Exist);
    end

% ======================================================================
% Function 2 : FM1: Full Information Forward Method
% ======================================================================
    %      Make agents use all variables of the the model as state variables 
    %      when forming expectations. 
    function [DETC,FCC,OmegaK,GammaK,FK]=fm1(A,B) % 
        tt=5; % Just in case the FCC holds, but with a large k such that termK=maxK
               % this step may be done more than once.
        for t1=1:tt
            if t1==1, H=0.1*ones(n,n); % 
            else, H=0.1*randn(n,n); %              
            end
                    [OmegaK,termK,FS_Exist]=fm1_OmegaK(A,B,H);  
                if FS_Exist==1,  break, end
            end
            
            [DETC,FCC,OmegaK,GammaK,FK]=fm_results(OmegaK,termK,FS_Exist);   
    end

    % Full Information Forward Algorithm for OmegaK to be controlled by fm1.m.
    
        function [OmegaK,termK,FS_Exist]=fm1_OmegaK(A,B,H)        
            FS_Exist=1;       
            % Select a well-defined Big matrices
         	rankBigB1=0; 
                while rankBigB1<2*n
                    BigB1_0=[eye(n) -eye(n);H*B eye(n)]; rankBigB1=rank(BigB1_0); 
                    if rankBigB1-2*n~=0, H=randn(n,n); BigB1_0=[eye(n) -eye(n);H*B eye(n)]; end
                end   
                BigB1=BigB1_0;
                BigA1=[zeros(n,2*n);A+H -H];
                BigB2=[B zeros(n,n);zeros(n,2*n)];
                BigA=BigB1\BigA1; BigB=BigB1\BigB2; N=2*n;
            % Initial 目_k at k=1
                BigOMEGAK{1,1}=BigB; 
            % 目_k for all k 
                k=0; FC_O=1;
                while FC_O>tolK
                    k=k+1;   
                        BigOMEGAK{k+1,1}=(eye(N)-BigA*BigOMEGAK{k,1})\BigB;
                    FC_O=max(abs(BigOMEGAK{k+1,1}(:)-BigOMEGAK{k,1}(:)));  
                    
                    % Checking whether FS is finite
                    if k==maxK, FS_Exist=0; break; end        
                    if max(abs(BigOMEGAK{k,1}(:)))>10e10, FS_Exist=0; break; end
                end 
            termK=k;   
            OmegaK=BigOMEGAK{end,1}(1:n,1:n);    
        end

% ======================================================================
% Function 3 : Modified Forward Method
% ======================================================================
    function [DETC,FCC,OmegaK,GammaK,FK]=fm(A,B) % 
    % Let FM0 (FM1) be Original (Full info) FM.
    % (1) Try FM0 first.
        [DETC,FCC,OmegaK,GammaK,FK]=fm0(A,B);                  
        if isnan(DETC(3)), [DETC,FCC,OmegaK,GammaK,FK]=fm1(A,B); 
    % (1-1) If DETC(3)<1 under FM0, Done.
    % (1-2) If DETC(3)>=1 or NaN, Apply FM1
        elseif DETC(3)>=1, [DETC_1,FCC_1,OmegaK_1,GammaK_1,FK_1]=fm1(A,B);
    % (1-2-1) If DETC(3)<1 under FM1, Done.
    % (1-2-2) If DETC(3)>=1 or NaN, Report Results that has smaller DETC(3)
    %         under FM0 and FM1.
            if DETC_1(3)<DETC(3),DETC=DETC_1; FCC=FCC_1; OmegaK=OmegaK_1; GammaK=GammaK_1; FK=FK_1;
            end
        end
    end

% ======================================================================
% Function 4 : Common Outputs from the computed OmegaK
% ======================================================================
    function [DETC,FCC,OmegaK,GammaK,FK]=fm_results(OmegaK,termK,FS_Exist)        
        if FS_Exist==1       
            Xi=eye(n)-A*OmegaK; 
            FK=Xi\A; 
            RtkFK=kron(R',FK);

            DETC1=max(abs(eig(OmegaK))); 
            DETC2=max(abs(eig(FK)));
            MODWD=DETC1*DETC2;
            FCC1=termK;
            FCC2=max(abs(eig(RtkFK)));  % Spectral Radius of kron(R',FK)        
            DETC=[DETC1 DETC2 MODWD]; %Determinacy Result+Condition for well-defined MOD solution
            FCC=[FCC1 FCC2]; 
            InvXiC=Xi\C;
            vGamma=(eye(n*m)-RtkFK)\InvXiC(:);
            GammaK=reshape(vGamma,n,m);
        elseif FS_Exist==0
            DETC=[NaN NaN NaN]; FCC=[NaN NaN]; 
            OmegaK=NaN*ones(n,n); GammaK=NaN*ones(n,m); FK=NaN*ones(n,n);
        end
    end
% ======================================================================
% end of Nested Function
% ======================================================================
         
        
end %% end of main function