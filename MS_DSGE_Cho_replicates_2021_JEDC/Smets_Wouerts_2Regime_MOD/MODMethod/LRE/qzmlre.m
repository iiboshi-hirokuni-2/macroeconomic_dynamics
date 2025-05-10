function [DETCMOD,OmegaMOD,GammaMOD,FMOD,Geig,varargout]=qzmlre(A,B,C,R)
%        [DETCMOD,OmegaMOD,GammaMOD,FMOD,Geig,gvAll,OmegaAll,GammaAll,FAll]=qzmlre(A,B,C,R) 
%% MAJOR REVISION of the previous code, written by Seonghoon Cho, August 8, 2018
% This is a matlab code implementing the MOD method for LRE model   
% using the eigensystem. See Cho(2019,WP).
%   (1) Compute the MOD solution using the eigensystem,
%   (2) Identify determinacy/indeterminacy/no stable solution (DET/INDET/NSS) 
%  Note:  This is comparable to gensys. Both use the same eigensystem.
%         This code modifies and incorporate the matlab codes written by Sims:
%         "reorder.m" and "qzswitch.m".
%  Can be downloaded from http://web.yonsei.ac.kr/sc719/ 
%% References :
%   Cho(2019,WP) "Determinacy and Classification of Markov-Switching Rational Expectations Models"
%% Usage: Specify the input arguments and run it
%       A (required) 
%       B,C,R and other optional arguments (optional).
%   The code will produce 
%       DETCMOD : Determinacy information+condition for determinacy-admissible model 
%       OmegaMOD,GammaMOD,FMOD: MOD Solution
%       Other output: 
%                    Geig  :Generalized Eigenvalues
%                    gvAll :Index of Geigs for all solutions
%                    OmegaAll,GammaAll,FAll: All solutions
%   This is a stand-alone code for implementing the MOD method.
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
%       (1) OmegaMOD = min r(Omega) for all Omega where r() is the spectral radius 
%           FMOD is the corresponding F  
%       (2) The model is Determinacy-Admissible iff 
%           OmegaMOD is real-valued and r(OmegaMOD)*r(FMOD)<1
%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       [1-1] Required 
%           A : n by n matrix (forward looking, may well be singular)
%       [1-2] Optional 
%           B : (Optional) n by n matrix (backwardlooking, may well be singular)
%                           Default : B=zeros(n,n);
%           C : (Optional) n by m coefficient matrix of the m by 1 vector of
%                exogenous variables, z(t). (Default: C=eye(n))
%           R : (Optional) VAR coefficient matrix of z(t-1) (Default:R=zeros(m,m))
%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DETCMOD :    [DETC1 DETC2 DETC1*DETC2]
%                       = [r(OmegaMOD) r(FMOD) r(OmegaMOD)*r(FMOD)] where r() is spectral radius.    
%           OmegaMOD :  Omega (MOD solution)
%           GammaMOD :  Gamma uniquely associated with OmegaMOD
%           FMOD :      F uniquely associated with OmegaMOD
%           Geig :      2n Generalized Eigenvalues 
%           gvAll :     N by n where N is number of all solutions
%                       each row is the set of n index out of 2n G-eigs.
%           OmegaAll :  N by 1 cell array: j-th element is j-th Omega  
%           GammaAll :  N by 1 cell array: j-th element is j-th Gamma
%           FAll :      N by 1 cell array: j-th element is j-th F        
%  
%      **** IMPORTANT NOTE. Computing all solutions can take a long
%           time because the total # of solutions can be very large 
%           reaching 184756 for n=10 and rank(B)=10.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Classification by the MOD Method
%                  Determinacy-Admissible                Determinacy-Inadmissible
%                  OmegaMOD is real-valued and           OmegaMOD is complex-valued or
%                  r(OmegaMOD)*r(FMOD)<1                 r(OmegaMOD)*r(FMOD)>=1 
%          DET     r(OmegaMOD)<1, r(FMOD)<=1             Impossible    
%          INDET   r(FMOD)>1                             r(OmegaMOD)<1     
%          NSS     r(OmegaMOD)>=1                        r(OmegaMOD)>=1  
%% INTERPRETATION of Results by the Eigensystem
% There are 2n generalized eigenvalues(G-eigs) associated with
% A and B where AA=[A 0;0 I] and BB=[I -B;I 0] are 2n by 2n matrices.
% G-eigs=[xi(1),xi(2),..,xi(n),xi(n+1),..,xi(2n)]' in ascending order in size.
% Omega is associated with a particular selection of n out of 2n G-eigs.
% Thus maximum # of Omega is (2n)!/[n!*n!]. 
% Definition
%   Omega_nSG = Omega associated with n smallest Geigs xi(1),xi(2),..,xi(n).
%   Determinacy-Admissible iff Omega_nSG exists AND |xi(n)|<|xi(n+1)| 
%% Classification by the Eigensystem
%                  Determinacy-Admissible                  Determinacy-Inadmissible
%                Omega_nSG exists AND |xi(n)|<|xi(n+1)|    Omega_nSG not exist or |xi(n)|=|xi(n+1)|
%          DET     |xi(n)|<  1 <=  |xi(n+1)|               Impossible    
%          INDET   |xi(n+1)| < 1                           |xi(n+i)|<1,   
%          NSS     1<= |xi(n)|                             1<= |xi(n+i)|
%% Summary: The following are equivalent:
%      1) MOD method
%      2) Root-Counting AND Existence of Omega_nSG
%      3) gensys algorithm  
%   Proof from a) r(OmegaMOD)= |xi(n)| and r(FMOD)=1/|xi(n+1)| if OmegaMOD =Omega_nSG
%              b) r(OmegaMOD)>=|xi(n)| and r(FMOD)>=1/xi(n)    if OmegaMOD~=Omega_nSG
%% Procedure of the Code
%  Check the existence of Omega_nSG.
%       If Omega_nSG exists, set OmegaMOD = Omega_nSG
%       If Omega_nSG does not exist, identify OmegaMOD by selecting next n smallest
%             Geigs.
% ======================================================================

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

%% [2] QZ Method   
%% [2-1] Apply Schur Decomposition and Arrange Geigs in increasing order.
        AA=[A zeros(n,n);zeros(n,n) eye(n)];
        BB=[eye(n) -B;eye(n) zeros(n,n)];
        
        [S0,T0,Q0,Z0] = qz(AA,BB);          % Initial Schur Decomposition (Not ordered)
        [S1,T1,Q1,Z1]=STQZ(S0,T0,Q0,Z0);    % Schur Decomposition in ascending order

        % Computing Generalized eigenvalues, Finite Geigs, Stable Geigs
            Geig=(diag(T0)./diag(S0))';
            for i=1:2*n
                if abs(Geig(i))>1e8, Geig(i)=sign(Geig(i))*Inf; end
                if abs(imag(Geig(i)))<1e-8, Geig(i)=real(Geig(i)); end
                if abs(S0(i,i))< 1e-8 && abs(T0(i,i))< 1e-8
                    error('*** Rank condition is violated (Coincident zeros): Model is not well-defined')
                end
            end 
            [~,indeig]=sort(abs(Geig));
            Geig=Geig(indeig);
            nf=sum(isfinite(Geig)); % # of finite Geigs  
%% [2-2] Computing MOD Solution = argmin sr(Omega) 
        % (1) Check the existence of the Omega_nSG=Omega(xi(1),xi(2),..,xi(n)).
                [S,T,Q,Z] = STQZgv(S1,T1,Q1,Z1,(1:1:n)',n); 
                [OmegaMOD,FMOD,GammaMOD,WD]=sol_gv(Z,n,m,A,C,R);  
                % If WD=1, then Omega_nSG exists and OmegaMOD=Omega_nSG. Done.
                
        % (2)   If ~Omega_nSG, do the following.
            %   - Identify maximum # of possible solutions
            %   - Find OmegaMOD such that min(r(Omega)).
            %     gvM(wdj,:) shows the order of Geigs of OmegaMOD.

                if WD==0
                    gvM=sortrows(nchoosek(1:1:nf,n),n); % Computing this can take time if n>=10.
                    ns=size(gvM,1); % maximum # of possible combinations of Geigs
                    wdj=0; %WD=0;
                    while WD<1 
                        wdj=wdj+1;
                        [S,T,Q,Z] = STQZgv(S1,T1,Q1,Z1,gvM(wdj,:)',n);
                        [OmegaMOD,FMOD,GammaMOD,WD]=sol_gv(Z,n,m,A,C,R);    
                        
                        if wdj==ns, break, end
                    end
                end    
               
        % (3) Checking DET/INDET/NSS 
                if isempty(OmegaMOD) || isempty(FMOD) 
                    DETC1=NaN; DETC2=NaN; % Won't happen
                else, DETC1=max(abs(eig(OmegaMOD))); DETC2=max(abs(eig(FMOD)));
                end

                DETCMOD=[DETC1 DETC2 DETC1*DETC2];            
            % Determinacy-Admissible if OmegaMOD is real-valued and DETC1*DETC2<1.    
           
%% [3] Other Output: All solutions 
        if nargout>5
                gvM=sortrows(nchoosek(1:1:nf,n),n); % Computing this can take time if n>=10.
                ns=size(gvM,1); % maximum # of possible combinations of Geigs
                kk=1; gvAll=[];
                for nsi=1:ns
                    [S,T,Q,Z] = STQZgv(S1,T1,Q1,Z1,gvM(nsi,:)',n);
                    [Omegaj,Fj,Gammaj,WDj]=sol_gv(Z,n,m,A,C,R);
                    if WDj==1
                        OmegaAll{kk,1}=Omegaj; FAll{kk,1}=Fj; GammaAll{kk,1}=Gammaj;
                        gvAll=[gvAll;gvM(nsi,:)];
                        kk=kk+1;
                    end              
                end     
                
            if ~exist('OmegaAll'), gvAll=[];  OmegaAll=[];  GammaAll=[]; FAll=[]; end  
               % Won't happen
               
            varargout(1)={gvAll}; varargout(2)={OmegaAll};
            varargout(3)={GammaAll}; varargout(4)={FAll};  
        end
        
%% End of the Main Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nested Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% ======================================================================
% Function 1 : Sorted Schur Decomposition
% ======================================================================
    function [S,T,Q,Z]=STQZ(S0,T0,Q0,Z0)
            nn = size(S0,1);  si = 1;
            while si<=nn-1
                if abs(T0(si,si)*S0(si+1,si+1))>abs(S0(si,si)*T0(si+1,si+1))   
                        % ====== Beginning of "qzswitch" block ==============
                        a = S0(si,si); d = T0(si,si); b = S0(si,si+1); e = T0(si,si+1);
                        c = S0(si+1,si+1); f = T0(si+1,si+1); 
                        wz = [c*e-f*b, (c*d-f*a)'];
                        xy = [(b*d-e*a)', (c*d-f*a)'];
                        sn = sqrt(wz*wz');
                        sm = sqrt(xy*xy');

                        if sn == 0, return         
                        else
                            wz = sn\wz;
                            xy = sm\xy;
                            wz = [wz; -wz(2)', wz(1)'];
                            xy = [xy;-xy(2)', xy(1)'];
                            S0(si:si+1,:) = xy*S0(si:si+1,:);
                            T0(si:si+1,:) = xy*T0(si:si+1,:);
                            S0(:,si:si+1) = S0(:,si:si+1)*wz;
                            T0(:,si:si+1) = T0(:,si:si+1)*wz;
                            Z0(:,si:si+1) = Z0(:,si:si+1)*wz;
                            Q0(si:si+1,:) = xy*Q0(si:si+1,:);
                        end
                        % ====== End of "qzswitch" block ==================
                    if ~(si==1);si = si-2;end
                end
                si=si+1;
            end
            S=S0;T=T0;Q=Q0;Z=Z0;
    end % end of function STQZ

% ======================================================================
% Function 2 : Schur Decomposition with a given gv 
% ======================================================================
    function [S,T,Q,Z]=STQZgv(S,T,Q,Z,gv,n)
        if isequal((1:1:n)',gv)==0
            if size(gv,1)==1, gv=gv'; end
            I=setdiff((1:1:n)',gv); % I=[i1,i2,...ik];
            J=setdiff(gv,(1:1:n)'); % J=[j1,j2,...jk];
            K=size(J,1);
               
            % Switch (ik,j1), (ik-1,k2)...(i1,jk)
            for k=1:K  
                ii=I(K+1-k); jj=J(k);
                
                for ki=ii:n-1
                    % ====== Beginning of "qzswitch" block ==============
                        a = S(ki,ki); d = T(ki,ki); b = S(ki,ki+1); e = T(ki,ki+1);
                        c = S(ki+1,ki+1); f = T(ki+1,ki+1); 
                        wz = [c*e-f*b, (c*d-f*a)'];
                        xy = [(b*d-e*a)', (c*d-f*a)'];
                        sn = sqrt(wz*wz');
                        sm = sqrt(xy*xy');

                        if sn ==0 , return         
                        else
                            wz = sn\wz;
                            xy = sm\xy;
                            wz = [wz; -wz(2)', wz(1)'];
                            xy = [xy;-xy(2)', xy(1)'];
                            S(ki:ki+1,:) = xy*S(ki:ki+1,:);
                            T(ki:ki+1,:) = xy*T(ki:ki+1,:);
                            S(:,ki:ki+1) = S(:,ki:ki+1)*wz;
                            T(:,ki:ki+1) = T(:,ki:ki+1)*wz;
                            Z(:,ki:ki+1) = Z(:,ki:ki+1)*wz;
                            Q(ki:ki+1,:) = xy*Q(ki:ki+1,:);
                        end
                    % ====== End of "qzswitch" block ==================          
                end

                for kj=1:jj-n
                    % ====== Beginning of "qzswitch" block ==============
                        a = S(jj-kj,jj-kj); d = T(jj-kj,jj-kj); b = S(jj-kj,jj-kj+1); e = T(jj-kj,jj-kj+1);
                        c = S(jj-kj+1,jj-kj+1); f = T(jj-kj+1,jj-kj+1); 
                        wz = [c*e-f*b, (c*d-f*a)'];
                        xy = [(b*d-e*a)', (c*d-f*a)'];
                        sn = sqrt(wz*wz');
                        sm = sqrt(xy*xy');

                        if sn ==0 , return         
                        else
                            wz = sn\wz;
                            xy = sm\xy;
                            wz = [wz; -wz(2)', wz(1)'];
                            xy = [xy;-xy(2)', xy(1)'];
                            S(jj-kj:jj-kj+1,:) = xy*S(jj-kj:jj-kj+1,:);
                            T(jj-kj:jj-kj+1,:) = xy*T(jj-kj:jj-kj+1,:);
                            S(:,jj-kj:jj-kj+1) = S(:,jj-kj:jj-kj+1)*wz;
                            T(:,jj-kj:jj-kj+1) = T(:,jj-kj:jj-kj+1)*wz;
                            Z(:,jj-kj:jj-kj+1) = Z(:,jj-kj:jj-kj+1)*wz;
                            Q(jj-kj:jj-kj+1,:) = xy*Q(jj-kj:jj-kj+1,:);
                        end
                    % ====== End of "qzswitch" block ==================
                end
          
            end % end of k=1:K; 
        end
    end % end of function STQZgv

% ======================================================================
% Function 3 : Solve for Omega, F and Gamma associated with given gv
% ======================================================================
    function [Omega,F,Gamma,WD]=sol_gv(Z,n,m,A,C,R)
        % [3-1] Defining Omega assiciated with gv a la McCallum
            H=inv(Z); H21=H(n+1:2*n,1:n);H22=H(n+1:2*n,n+1:2*n);
            
            if rank(H21)==n % Omega to be well-defined 
                % Compute Omega
                    Omega=(-inv(H21)*H22);
                        for si=1:n, for sj=1:n
                            if abs(imag(Omega(si,sj)))<1e-8, Omega(si,sj)=real(Omega(si,sj)); end
                        end, end
                % Compute F and Gamma
                    Xi=eye(n)-A*Omega;    
                    if rank(Xi)==n % F to be well-defined 
                        F=Xi\A;
                        RtkFK=kron(R',F);
                        InvXiC=Xi\C;
                        vGamma=(eye(n*m)-RtkFK)\InvXiC(:);
                        Gamma=reshape(vGamma,n,m);
                    else, F=[]; Gamma=[]; Omega=[]; 
                    end
                    
            else, Omega=[]; F=[]; Gamma=[];
            end
        % [3-2] Identifying whether the solution is well-defined.
                if isempty(Omega) || isempty(F), WD=0;  
                else, WD=1; 
                end
      
    end  % end of function sol_gv    
% ======================================================================
% end of Nested Function
% ======================================================================

end % End of Main Function
