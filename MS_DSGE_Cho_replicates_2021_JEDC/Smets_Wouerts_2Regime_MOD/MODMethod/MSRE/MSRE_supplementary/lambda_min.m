function [Lmin,tau2,xi2,V,Phi,D,u,Q]=lambda_min(P,F)
 
%% [1] Identify S, n from F. Construct Phi_kron(F,F) and compute xi2:
    S=size(P,1);
    n=size(F{1,1},1);
    N=n^2*S;

    [PhiK_F_F,eigPhiK_F_F,xi2]=PHIK(P,F,F);
    
%% [2] Transformation of the model with Fh and Ph such that 
%%     Ph = ones(S,S),  Fh=F^hat and r(Phi_kron(Fh,Fh))=1.
    Fh=cell(S,S); for i=1:S, for j=1:S, Fh{i,j}=F{i,j}*sqrt(P(i,j)/xi2); end, end
    Ph=ones(S,S);    
    
    [PhiK_Fh_Fh,eigPhiK_Fh_Fh,xi2_Fh]=PHIK(Ph,Fh,Fh);
            

%% [3] Finding i_all, the set of the indices of the eigenvector u such that 
%%           PhiK_Fh_Fh * u = u
    [Uh,Dh]=eig(PhiK_Fh_Fh); [~,Idxh]=max(abs(diag(Dh))); 
    u=Uh(:,Idxh);
       
%% [4] Define Q(s)=reshape(u(s)) where u=[u(1);...;u(S)] 
%%     Define n by n VA(s) and DA(s) such that Q(s)VA(s)=VA(s)DA(s) 
%%     using Schur decomposition.
%%     Select a set of non-zero eigenvalues of DA such that
%%     Q(s)V(s)=V(s)D(s) where V(s) is n by ks and D(s) is ks by ks.
%%       0<=ks<=n  
     
    us=cell(S,1); % Subvector of u at each regime s
    Q=cell(S,1);  V=cell(S,1); D=cell(S,1); rD=cell(S,1); 
        for s=1:S
            us{s,1}=u(n^2*(s-1)+1:n^2*s,1);
            
            Q{s,1}=reshape(us{s,1},n,n); 
            [Vtmp,Dtmp]=schur(Q{s,1});
            for i=1:n, for j=1:n
                    
                    if abs(Dtmp(i,j))<10^(-6), Dtmp(i,j)=0; end % Detecting 0 by removing computational errors. 
            end, end
            DE=ordeig(Dtmp);
            
            [Vtmp,Dtmp]=ordschur(Vtmp,Dtmp,abs(DE)>0.0000001); % Selecting eigenvectors associated with non-zero eigenvalues.
            rD{s,1}=rank(Dtmp);
            V{s,1}=Vtmp(:,1:rD{s,1});      
            D{s,1}=Dtmp(1:rD{s,1},1:rD{s,1});          
        end
        
%% [5] Defining D(s) and V(s) 
        for i=1:S, for j=1:S
                Phi{i,j}=D{j,1}*V{j,1}'*F{i,j}'*V{i,1}*D{i,1}^(-1)/xi2;
                Lmin{i,j}=V{j,1}*Phi{i,j}*V{i,1}';
            end
        end

        [BarPhiK_L_L,eigBarPhiK_L_L,tau2]=BarPHIK(P,Lmin,Lmin);
        
if abs(xi2*tau2-1)>0.00005
    disp('===============================================================')
    disp('xi2*tau2 is') 
    disp(xi2*tau2)
    disp('This may be because the choice of D and V is not correct.')
    disp('Lower precision to choose alternative D and V and run again.')
    disp('Or because the matlab code eig/schur may not be') 
    disp('accurate when A=BarPhiK_L_L is ill-conditioned. In this case')
    disp('eig(A) may not be the same as eig(A^(T)). ')
end


%% Nested Functions
% ======================================================================
function [PhiK_XY,eigPhiK_XY,mai]=PHIK(P,X0,Y0,tp1,tp2)

S=size(P,1); 

if iscell(X0)==0, n=size(X0,1); for i=1:S, for j=1:S, X{i,j}=X0; end, end
else, n=size(X0{1,1},1); X=X0;
end

if iscell(Y0)==0, for i=1:S, for j=1:S, Y{i,j}=Y0; end, end
else, Y=Y0;
end

if size(X,2)==1, for i=1:S, for j=1:S, X{i,j}=X{i,1}; end, end, end
if size(Y,2)==1, for i=1:S, for j=1:S, Y{i,j}=Y{i,1}; end, end, end

if nargin>=4 && tp1==1
    for i=1:S, for j=1:S, X{i,j}=X{i,j}'; end, end
end

if nargin==5 && tp2==1
    for i=1:S, for j=1:S, Y{i,j}=Y{i,j}'; end, end
end

    PhiK_XY=[]; 
        for i=1:S, tmpr=[];
            for j=1:S, tmpr=[tmpr P(i,j)*kron(X{i,j},Y{i,j})]; end
            PhiK_XY=[PhiK_XY;tmpr];
        end
        tmp=eig(PhiK_XY)'; [~,idx]=sort(abs(tmp),'descend');
        eigPhiK_XY=tmp(1,idx);
    mai=max(abs(eig(PhiK_XY)));  
end

% ======================================================================
function [BarPhiK_XY,eigBarPhiK_XY,mai]=BarPHIK(P,X0,Y0,tp1,tp2)

S=size(P,1); 

if iscell(X0)==0, n=size(X0,1); for i=1:S, for j=1:S, X{i,j}=X0; end, end
else, n=size(X0{1,1},1); X=X0;
end

if iscell(Y0)==0, for i=1:S, for j=1:S, Y{i,j}=Y0; end, end
else, Y=Y0;
end

if size(X,2)==1, for i=1:S, for j=1:S, X{i,j}=X{i,1}; end, end, end
if size(Y,2)==1, for i=1:S, for j=1:S, Y{i,j}=Y{i,1}; end, end, end

if nargin>=4 && tp1==1
    for i=1:S, for j=1:S, X{i,j}=X{i,j}'; end, end
end

if nargin==5 && tp2==1
    for i=1:S, for j=1:S, Y{i,j}=Y{i,j}'; end, end
end

BarPhiK_XY=[]; 
        for i=1:S, tmpr=[];
            for j=1:S, tmpr=[tmpr;P(i,j)*kron(X{i,j},Y{i,j})]; end
            BarPhiK_XY=[BarPhiK_XY tmpr];
        end
    tmp=eig(BarPhiK_XY)'; [~,idx]=sort(abs(tmp),'descend');
    eigBarPhiK_XY=tmp(1,idx);
mai=max(abs(eig(BarPhiK_XY)));  


end


end