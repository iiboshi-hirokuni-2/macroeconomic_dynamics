% 
%    solving Markov Swicht DSGE 
%      by Farmer Waggoner and Zha (2009, 2011)
% 
% 

warning 'off'

% close all

disp( '  '  )
disp('=============================')
disp('           start MS DSGE model ')
disp('=============================')

%% setting 

 %  transtion Probs
  
  MS_transition_prob;

 % 
 neta =3;  % Numbers of  Jump variables

ns = size(PP,1);
nvars=neq ;
nexperrs =neta;
max_iterations_msv_msre = 100; 
smallval = 3.25; %1.0e-1;  %Convergence criterion: sqrt(machine epsilon).
  ini_scale = 7.5e-3 ;

disp(' ')
disp('Ergodic probability:')
fn_ergodp(PP)

%%
% SET MATRIX for FWZ method
i = 1; % regime 1
     Afwz_cell{i}=GAM01;
     Bfwz_cell{i}=GAM11;
     gPsifwz_cell{i}=PSI01;
     PPI(:,:,i)=PPI1;

% SET MATRIX for FWZ method
i = 2;  % regime 2 
    Afwz_cell{i}=GAM02;
     Bfwz_cell{i}=GAM12;
     gPsifwz_cell{i}=PSI02;
     PPI(:,:,i)=PPI2;

%---------------------------------------------------------------------------
%------------------------ $$ The new FWZ algorithm $$ ----------------------
%---------------------------------------------------------------------------

geneigvals = zeros(nvars,ns);

ngensys = size(Afwz_cell{1},1);

rng(100)  % random number  seed()

[F1fwz_cell, F2fwz_cell, G1fwz_cell, G2fwz_cell, Vfwz_cell, flag_err] ...
                                       = fwz_msv_msre(PP', Afwz_cell, Bfwz_cell, gPsifwz_cell, nexperrs,...
                                                                    ini_scale*randn(ns*nexperrs*(ngensys-nexperrs),1),...
                                                                    max_iterations_msv_msre,smallval);

 %---------------------------------------------------------------------------------------                                     
 %---------------------------------------------------------------------------------------                                          

 if (flag_err<0)
   disp('  ')
   disp(' ')
   disp('------------- There is no MSV equilibrium! ------------')
   disp(' ')
else
   disp('  ')
   disp(' ')
   disp('------------- FWZ solution ------------')
   disp('Solution G1(j) for x_t = G1(j) x_{t-1} + G2(j) \epsilon_t,')
   disp('  where x_t = [x_t, E_t x_{t+1}, pi_t, E_t pi_{t+1}, R_t, z_{D,t}, z_{S,t}]'' for both regimes:')
   Gfwz_cell = cell(ns,1);
   Impactfwz_cell = cell(ns,1);
   for si=1:ns
      Gfwz_cell{si} = Vfwz_cell{si} * F1fwz_cell{si};
       Gfwz_cell{si}
   end   
   disp('Impact matrix G2(j) on \epsilon_t for both regimes:')
   for si=1:ns
      Impactfwz_cell{si} = Vfwz_cell{si} * G1fwz_cell{si};
      Impactfwz_cell{si}
   end
   %====== Checking stationarity of the solution. =======
   %   Checking the stationarity of regime-switching models. 
   %     If the rows of P sum to one and Gamma(i)=V(i)*A(i) is n x n, then the solution is MSS stable if and only if the eigenvalues of
   %            kron(P',eye(n*n)) *diag(kron(Gamma(i),Gamma(i)))
   %     are all inside the unit circle.
   nfwz = size(Afwz_cell{1},1);
   stackblkdiag = zeros(ns*nfwz^2);
   for si=1:ns
      stackblkdiag(((si-1)*nfwz^2+1):(si*nfwz^2),((si-1)*nfwz^2+1):(si*nfwz^2)) = kron(Gfwz_cell{si},Gfwz_cell{si});
   end 
   stackmat = kron(PP,eye(nfwz*nfwz)) * stackblkdiag; 
   disp('------------ ')
   disp('If the following (max) absolute value of eigenvalues is less than 1.0, then stationary.')
   max(abs(eig(stackmat)))
   eig(stackmat);
end                                     

ngensys = size(Afwz_cell{1},1);
nsh = size(gPsifwz_cell{1},2);

    T1=zeros(ngensys,ngensys,h);
    T0=zeros(ngensys,nsh,h);
    for si=1:h
        T1(:,:,si) = Vfwz_cell{si} * F1fwz_cell{si};
        T0(:,:,si) = Vfwz_cell{si} * G1fwz_cell{si};
    end   


pause(0.5)


%%
% Plot_IRF_MS_R2;
    