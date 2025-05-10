function [DETCMOD,OmegaMOD,FMOD,DETC_All,AllOmegas,DETC_OmjtkF1]=gbmsre(P,A,B) 
%% This code implements the MOD method of Cho(2019) using the Groebner basis(GB) approach.
%   - The code call "Singular" program in actual computation of the solutions, thus 
%   - Users need to do the following.
%% Usage 
%   (1) Install Singular
%   (2) Write your n-dimensional MSRE model with S regimes in terms of P, A and B and run it.
%      : Computation may fail even for a very small scale model. Refer to Testing_GB_computation_time.m.
%   (3) Follow the instruction displayed in matlab command.
%   (4) Please read carefully the instruction below.
%% Acronyms
%   FM /GB : Forward Method / Groebner Basis Approach
%   DA / DIA : Determinacy-Adimissible / Determinacy-Inadmissible
%   DET/INDET/NSS : Determinate / Indeterminate / No Stable Solution
%   MOD : MOD solution
%% ================================================================
%% Preparations : Please read carefully the following. 
%% (1) Install Singular from the following website at default location in C drive. 
%%     https://www.singular.uni-kl.de/index.php/singular-download.html
%%     Open Cygwin64 Terminal and type 'singular' in that command window.
%    	- "C:\cygwin64\home\Computername" is the folder in which this matlab code
%          writes the Singular procedure code, Singular reads the code and computes
%          and store the solutions as ascii file, and this matlab code load the file
%          and store the results in the current folder you work at. 
%      	-  Your Computer name will be automatically detected, so you don't need to do any thing.
%      	-  See Section [0] below if error message pops up.
%% (2) Wite your own matlab code describing the model in terms of P, A, B. 
%      	- See Section [1]
%% ================================================================
%% Understanding How this code works.
%% (1) Generating a code to be read by Singular.
%      - This code writes a Singular execution code "GB_singular.ascii",
%        which transforms the solution restrictions implied by your model into a set of polynomials.
%        Then this code pauses temporarilly.
%% (2) Ask Singular to compute the all of the solutions using GB.
%%     To do so, copy and paste " execute(read("GB_singular.ascii")); " 
%%     in Singular command prompt and press Enter.
%   	- The sinular code will compute all of the solutions for Omega
%       - Store them as Sol_Rpart.txt and Sol_Cpart.txt,
%       - which are the real and imaginary parts of the solutions.
%       - *** IMPORTANT: The GB algorithm may run forever if the model is 
%         more than 3 dimensional with more than two lagged variables.
%       - If you see a prompt ">" in the command line of Singular, computation is done. 
%% (3) Go back to matlab command window and press any key.
%       - This code will load Sol_partR.txt and Sol_partC.txt and convert them into
%         the output of this file. The results will be displayed. 
%% ================================================================
%% STRUCTURE OF THIS CODE
%%      [Section 0] : Adding your sygwin directory to filepath of matlab.
%%      [Section 1] : Handling Inputs
%%      [Section 2] : Writing a Singular code 
%%      [Section 3] : Ask Singular to solve the problem 
%%                    Your action needed: Follow instruction in matlab command window
%%      [Section 4] Transform the Singular results into the output format in this code.
%%      [Section 5] Display Results
%%      [Section 6] Auxiliary Functions
%% 
%% ================================================================
%% NOTES 
%% (1) This code is a preliminary version of implementing the MOD method
%%     using matlab and Singular at the same time.
%%     Will be updated if necessary. Comments will be greatly appreciated.
%% (2) Written by Seonghoon Cho, Mar 14, 2019. 
%% ==============================================================
%% INPUT
%       P,A,B [Required]: Same input as those used in fmmsre.m : Refer to the code.
%% OUTPUT 
%       DETCMOD : The first row of DETC_All : Information for MOD solution
%       OmegaMOD, FMOD : MOD solution and the corresponding F : Same format as OmegaK, FK in fmmsre.m 
%       DETC_All  Ns by 5 matrix :  for each row solution Omegaj and Fj, 
%               DETC1=r(BarPsi_kron(Omegaj,Omegaj)) : 
%               DETC2=r(Psi_kron(Fj,Fj)).
%               DETC3=DETC1*DETC2
%               DETC4=r(Psi_kron(Omegaj',Fj));
%               DETC5=1 if the solution is "real-valued" and 0 otherwise.
%       AllOmegas: ns by S cell array where ns is the total number of solutions
%               AllOmegas(j,s) is the j-th solution Omega(s(t)=s). 
%       DETC_OmjtkF1=r(Psi_kron(Omegaj',F1)). This shows an additional result
%               for Proposition 4 to show that there 
%               exists Omegaj such that r(Psi_kron(Omegaj',F1))=1.
%   Note : DETC_All and AllOmegas are arranged in an increasing order of DETC1.
       
%% [Section 0] : Checking your home directory of Singular.
        listing=dir("C:\cygwin*"); 
            if ~isempty(listing) && listing.isdir==1, dirname=listing.name; 
            else, disp('Error: make sure you installed Singular in C:\cygwin*  '), return
            end 
        tmp_dir=strcat('C:\',  dirname,  '\home\'); 
        MyFolderInfo = dir(tmp_dir); MyFolderInfo.name;
        home_dir=strcat(tmp_dir, MyFolderInfo(3).name); 
        addpath(home_dir);

        outputfile="GB_singular.ascii"; 
         
%% [Section 1] Handling Inputs
    S=size(P,1);        % number of regimes
    if ~iscell(A), aa=A; clear A; A{1,1}=aa; end
    if ~iscell(B), bb=B; clear B; B{1,1}=bb; end
    n=size(A{1,1},1);   % Dimension of the model
    N=n^2*S;            % Total number of unknowns in Omega(s(t))
    if size(A,2)==1, AA=cell(S,S); 
        for i=1:S, for j=1:S, AA{i,j}=A{i,1}; end, end, A=AA; 
    end
    
%% [Section 2] Writing a Singular code implementing Groebner Basis Approach
%% [2-1] Generating Polynomials associated with the model in Singular form
    % [2-1-1] : Define N polynomials (xv_poly) in unrestricted N unknowns.
        X=cell(S,1);   % Define X=Omega, the ij-th entry of Omega(s) as x_s_ij
            for Si=1:S
            X_Si_name=strcat("x_", string(Si), "_%d%d");
            X{Si}=sym(X_Si_name, [n n]);
            end
        xv=X{1}(:); % vectorize [Omega(s(t)=1);...Omega(s(t)=S)]
                     % This is N by 1 unknowns
            for Si=2:S, xv=[xv;X{Si}(:)]; end      
        
        X_poly=cell(S,1);  % X_poly(s(t))=E[A(s(t),s(t+1))Omega(s(t+1))]Omega(s(t)-Omega(s(t)+B(s(t))=0
        for Si=1:S
            tmp_sum_j=P(Si,1)*A{Si,1}*X{1};
            for Sj=2:S, tmp_sum_j=tmp_sum_j+P(Si,Sj)*A{Si,Sj}*X{Sj}; end
            X_poly_Si=tmp_sum_j*X{Si}-X{Si}+B{Si,1};
            X_poly{Si}=X_poly_Si;
        end

        xv_poly=X_poly{1}(:); % vectorize [X_poly(s(t)=1);...X_poly=S)]
                      % This is N by 1 multivariate polynomials
        for Si=2:S, xv_poly=[xv_poly;X_poly{Si}(:)];  end
        
    % [2-1-2] : Identify true N_nz state variables and N_z non-state
    %           variables out of N unknowns in Omega(s(t)).
        xv_z=xv;  
        % (Case1) Identify the i-th columns of zeros in B(s(t)) for all s(t).
        B_zero_idx=cell(S,1); 
        for Si=1:S, B_zero_idx{Si,1}=ones(n,n); end
        Btmp=abs(B{1,1}); 
        for Si=2:S, Btmp=Btmp+abs(B{Si,1}); end
        for j=1:n 
                if Btmp(:,j)==zeros(n,1) 
                    for Si=1:S, B_zero_idx{Si,1}(:,j)=zeros(n,1); end
                end
        end
        
        v_zero_idx=B_zero_idx{1,1}(:);
        for Si=2:S, v_zero_idx=[v_zero_idx;B_zero_idx{Si,1}(:)]; end
        i_z=[];  
        i_nz=[];% Index of true unknowns out of N.
        for Ni=1:N 
            if v_zero_idx(Ni,1)==0,xv_z(Ni)=0; i_z=[i_z;Ni]; 
            else, i_nz=[i_nz;Ni]; 
            end 
        end
        xv_z(i_z,1)=zeros(size(i_z,1),1);
        
   % [2-1-3] Redefine N_nz polynomials (xv_nz_poly) in N_nz unknowns.
        xv_z(i_z)=double(xv_z(i_z));
        X_z=cell(S,1);
        for Si=1:S, X_z{Si}=reshape(xv_z((Si-1)*n^2+1:Si*n^2,1),n,n); end
        
        X_z_poly=cell(S,1);
        Omega_stack=X_z{1}; for Si=2:S, Omega_stack=[Omega_stack;X_z{Si}]; end
        for Si=1:S 
            tmp_stack_j=P(Si,1)*A{Si,1}; 
            for Sj=2:S, tmp_stack_j=[tmp_stack_j P(Si,Sj)*A{Si,Sj}]; end
            tmp_stack_j=round(tmp_stack_j,5);
            X_z_poly_Si=tmp_stack_j*Omega_stack*X_z{Si}-X_z{Si}+round(B{Si,1},5);
            X_z_poly{Si}=X_z_poly_Si;
        end
        
        xv_z_poly=X_z_poly{1}(:);
        for Si=2:S, xv_z_poly=[xv_z_poly;X_z_poly{Si}(:)];  end
          
        xv_nz=xv_z(i_nz,1);
        xv_nz_poly=xv_z_poly(i_nz,1);
        N_nz=size(xv_nz_poly,1);
             
%% [2-2] Writing the code Implementing the GB method in Singular
        
        f = sym('f', [1 N_nz]); % total number of non-zero polynomials. 
        strR="Sol_partR.txt"; 
        strC="Sol_partC.txt";
        strRw=strcat("write("":w ",strR,""",SOLMR);");
        strCw=strcat("write("":w ",strC, """,SOLMC);");
        
        st{1,1}='timer=0; system("--ticks-per-sec",1000); int t=timer; // CPU time in milliseconds.\n';
        
        %st{1,1}='timer=1; int t=timer; // CPU time in seconds;\n';
        
        st{2,1}='LIB "solve.lib";\n';
        st{3,1}='ring rr=(real,i),xxx,dp;\n';
        st{4,1}='// Beginning of Model Specification\n';
            str_n=strjoin(["int n=" string(n) "; // model dimension\n"]);
        st{5,1}=char(str_n);
            str_S=strjoin(["int S=" string(S) "; // number of regimes\n"]);
        st{6,1}=char(str_S);
        st{7,1}='int N=n^2*S; // total number of unknowns\r\n';
            str_N_nz=strjoin(["int N_nz=" string(N_nz) "; // number of non-zero unknowns\n"]);
        st{8,1}=char(str_N_nz);
            ringdef=["ring r=0,("];
            for i=1:N_nz-1 
                ringdef=[ringdef string(xv_nz(i)) ","]; 
            end
            ringdef=[ringdef string(xv_nz(N_nz)) "),dp;\n"];
        st{9,1}=char(strjoin(ringdef));
            strex=[];
            for j=1:N_nz
                strex=[strex;strjoin(["poly", string(f(j)),"=", string(xv_nz_poly(j)) ";\n"])];
            end       
            n_strex=size(strex,1);
            for j=1:n_strex,
                st{9+j,1}=char(strex(j,1));
            end
            idealdef=["ideal fi="];
            for j=1:N_nz-1
                idealdef=[idealdef string(f(j)) ","]; 
            end
            idealdef=[idealdef string(f(N_nz)) ";\n"];
            Newn=9+n_strex;     
        st{Newn+1,1}=char(strjoin(idealdef));
        st{Newn+2,1}=' \n';
        st{Newn+3,1}= '// option(prot);\n';
        st{Newn+4,1}='ideal si=groebner(fi);\n';
        st{Newn+5,1}='int nc=ncols(si);\n';
        st{Newn+6,1}='int ns=vdim(si); // maxinum number of solutions\n';
        st{Newn+7,1}='setring r;\n';
        st{Newn+8,1}='def T=solve(si,6,1,"nodisplay"); // Solve for all solutions using Groebner Basis;\n';
        st{Newn+9,1}='setring T;\n';
        st{Newn+10,1}=' \n';
        st{Newn+11,1}='// Store Solutions into Matrices;\n';
        st{Newn+12,1}='matrix SOLMR[ns][N_nz];\n';
        st{Newn+13,1}='int ii,jj;\n';
        st{Newn+14,1}='for (ii=1;ii<=ns;ii++)\n';
        if N_nz>1
            st{Newn+15,1}='{  for (jj=1;jj<=N_nz;jj++) { SOLMR[ii,jj]=repart(SOL[1][1][ii][jj]); } };\n';
        elseif N_nz==1
            st{Newn+15,1}='{  SOLMR[ii,1]=repart(SOL[1][1][ii]); };\n';
        end
        st{Newn+16,1}=' \n';
        st{Newn+17,1}='matrix SOLMC[ns][N_nz];\n';
        st{Newn+18,1}='int ii,jj;\n';
        st{Newn+19,1}='for (ii=1;ii<=ns;ii++)\n';
        if N_nz>1
            st{Newn+20,1}='{  for (jj=1;jj<=N_nz;jj++) { SOLMC[ii,jj]=impart(SOL[1][1][ii][jj]); } };\n';
        elseif N_nz==1
            st{Newn+20,1}='{  SOLMC[ii,1]=impart(SOL[1][1][ii]); };\n';
        end
        st{Newn+21,1}=strjoin([strRw "\n"]);
        st{Newn+22,1}=strjoin([strCw "\n"]);
        st{Newn+23,1}=' // Go back to matlab command window and press any key.;\n';
        st{Newn+24,1}=' \n';

        outputpath=strcat(home_dir,"\",outputfile);
        fileID = fopen(outputpath,'w');
            for i=1:size(st,1)
                fprintf(fileID,st{i});
            end
        fclose(fileID);



%% [Section 3] Ask Singular to solve the problem     
        disp('   ')
        disp('*** IMPORTANT !!! : DO THE FOLLOWING. ***')    
        disp('(1) Open the Singular window: Run Cygwin64 Terminal and type "singular" ')
        disp('(2) Copy the following and paste it in the Singular command line and')
        disp('    Press enter key and wait until the code finish computation.')
            commandtext= strjoin(['        ',strcat('execute(read("',outputfile, '")); timer; ns;')]); 
            disp(commandtext)
        disp('(3) Wait until a number shows up, which means that computation is finished.')
        disp('    The numbers are the time elapsed in milliseconds, and the number of solutions.')
        disp('    Come back to matlab window and press any key.')

    pause % Wait until Singular finishes computation.

%% [Section 4] Transform the Singular results into the output format in this code.
%% [4-1] Retrieving the solutions from GB       
        tmp_solR=load(strR);
        tmp_solC=load(strC);
        tmp_sol=complex(tmp_solR,tmp_solC);
        TN_nz=size(tmp_sol,2);
        ns=TN_nz/N_nz; % # of Solutions
%% [4-2] Write the solutions into the original form of Omega(s(t))  
        solmt=reshape(tmp_sol,N_nz,ns).'; % ns by N_nz matrix.
        AllOmega=cell(ns,S);
        for nsi=1:ns
            vsol=zeros(1,n^2);
            vsol(1,i_nz)=solmt(nsi,:);
            if ~isempty(i_z),
                vsol(1,i_z)=double(xv_z(i_z,1))';
            end
            vSsol=reshape(vsol,N/S,S);
            for Si=1:S               
                AllOmega{nsi,Si}=reshape(vSsol(:,Si),n,n);
            end
        end
    
%% [4-3] Computing all solutions and DETC information

        [DETCMOD,OmegaMOD,FMOD,DETC_All,AllOmegas,DETC_OmjtkF1]=gb_results(P,A,AllOmega); 
        % Sorting the solutions in an increasing order of DETC(1).
        Sol_order=string((1:1:ns)');
        Sol_order(1)="MOD"; 
        if isequal(DETC_All(1,1),DETC_All(2,1)), Sol_order(2)="MOD"; end

        mdi=1; % Do not display if DETC(3)>>1000.
        for i=1:ns,
            if max(DETC_All(i,:))>=1000, mdi=i-1; break, 
            else, mdi=ns; 
            end 
        end

        D5str=string(DETC_All(:,5)); 
        for nsi=1:ns
            if DETC_All(nsi,5)==1, D5str(nsi,1)='real_valued'; 
            elseif DETC_All(nsi,5)==0, D5str(nsi,1)='complex_valued';
            end
        end
        
%% [Section 5] Display Results
        disp('   ')
        disp('=====================================================')  
        disp('===== RESULTS from Groebner Basis Approach ==========')  
        disp('=====================================================')  
        disp(strjoin(["    Total number of solutions = " ns ]))
        disp('    Table below shows DETC(1) through DETC(4) information and real-valuedness')  
        disp('    for solutions i=1:ns with DETC(3)<1000.')
        disp('    ')
        varNames = {'Solution_Order','DETC','Solution_is'};
        T2=table(Sol_order(1:mdi,1),DETC_All(1:mdi,1:4),D5str(1:mdi,1),'VariableNames',varNames);
        disp(T2)
        disp('    The MOD solution  Omega(s) is')
        for s=1:S, disp(AllOmegas{1,s}), end
        
        %    disp('    If the MOD solution is real-valued, and DETC1(Omega^{2})*DETC2(OmegaMOD)>1,')
        %    disp('    this is an example that uniqueness of stable MSV solution does not imply Determinacy.')
        %    disp('    This may be directly such an example. If not, set for this model ') 
        %    disp('     alpha=(DETC_All(1,2)/DETC_All(2,1))^(1/4),')
        %    disp('    Define A( ) = A( )/alpha; B(  )=alpha*B(  ) for all regimes, and run this again.')
        %    disp('    Refer to Proposition 4 of Cho(2019, WP)')

%% [Section 6] Auxiliary Functions
function [DETCMOD,OmegaMOD,FMOD,DETC_All,AllOmegas,DETC_OmjtkF1]=gb_results(P,A,AllOmega) 
[ns,S]=size(AllOmega);
n=size(AllOmega{1,1},1);
N=n^2*S;

DETC_All=zeros(ns,5); 
for nsi=1:ns
    Omega_nsi=AllOmega(nsi,:);
    F_nsi=Ffun(Omega_nsi,A);
    BarPsi_OmkOm_nsi=BarPsi_X1X1(Omega_nsi);
    Psi_FkF_nsi=Psi_XX(F_nsi);
    Psi_OmtkF_nsi=Psi_X1tY(Omega_nsi,F_nsi);
    DETC_All(nsi,1)=max(abs(eig(BarPsi_OmkOm_nsi)));
    DETC_All(nsi,2)=max(abs(eig(Psi_FkF_nsi)));
    DETC_All(nsi,3)=DETC_All(nsi,1)*DETC_All(nsi,2);
    DETC_All(nsi,4)=max(abs(eig(Psi_OmtkF_nsi)));
        tmpreal=zeros(1,S); 
        for Si=1:S, tmpreal(1,Si)=isreal(Omega_nsi{Si}); end
    if tmpreal==ones(1,S), DETC_All(nsi,5)=1; 
    else, DETC_All(nsi,5)=0; 
    end
end

[DETC_All,sort_idx]=sortrows(DETC_All,[1 3]);
AllOmegas=AllOmega(sort_idx,:);
DETCMOD=DETC_All(1,:);
OmegaMOD=AllOmegas(1,:)'; FMOD=Ffun(OmegaMOD,A);

% Checking whether there exists j such that r(Psi_kron(Omega_MOD',Fj))=1.
DETC_OmjtkF1=zeros(ns,1);
for nsi=1:ns
    Omega_nsi=AllOmegas(nsi,:);
    DETC_OmjtkF1(nsi,1)=max(abs(eig(Psi_X1tY(Omega_nsi,FMOD))));
end
    
    
%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Computing F from Omega
    function F_i=Ffun(Omega_i,A)
            for i=1:S
                EAOmega_i=zeros(n,n); 
                    for j=1:S               
                        EAOmega_i =EAOmega_i +P(i,j)*A{i,j}*Omega_i{j};
                    end
                    Xi_i{i,1}=eye(n)-EAOmega_i;
            end 
            for i=1:S, for j=1:S, F_i{i,j}=Xi_i{i,1}\A{i,j}; end, end 
    end
    %% BarPsi_kron(Omega,Omega)
    function BarPsi_XkX=BarPsi_X1X1(X)
        bdiagXX=zeros(n^2*S,n^2*S);
        for i=1:S
            bdiagXX(n^2*(i-1)+1:n^2*i,n^2*(i-1)+1:n^2*i)...
            = kron(X{i},X{i});
        end
        BarPsi_XkX= bdiagXX*kron(P',eye(n^2));  %
    end
    
    %% Psi_kron(F,F)
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
    %% Psi_kron(Omega',F)
    function Psi_X1tkY=Psi_X1tY(X1,Y)
        Psi_X1tkY=[]; 
        for i=1:S
            Psi_X1tkYrow=[]; 
            for j=1:S
                    Psi_X1tkYrow=[Psi_X1tkYrow P(i,j)*kron(X1{j}',Y{i,j})];
            end
            Psi_X1tkY=[Psi_X1tkY;Psi_X1tkYrow];  
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end