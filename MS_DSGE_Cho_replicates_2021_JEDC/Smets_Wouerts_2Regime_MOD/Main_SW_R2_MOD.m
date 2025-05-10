% disp( '  '  )
% disp('============================================')
% disp('           start  Smets & Wouters (AER 2007) DSGE model ')
% disp(          ' Solution method ; MOD version ')
% disp('============================================')

clc

addpath(genpath('./MODMethod'))
addpath('./Func/')

clear all
close all

disp( '  '  )
disp('============================================')
disp('           start  Smets & Wouters (AER 2007) DSGE model ')
disp(' Solution method ; MOD version ')
disp('============================================')

%% Load  Smets & Wouters (AER 2007) DSGE model
SW_model;

%%
B1j_R1 = -1*jacobian(System_of_Eq, Y_t);
A1j_R1 = jacobian(System_of_Eq, Y_f);
B2j_R1  = jacobian(System_of_Eq, Y_b);
C1j_R1 = jacobian(System_of_Eq, Epsilon_t);

B1j_R2 = -1*jacobian(System_of_Eq, Y_t);
A1j_R2 = jacobian(System_of_Eq, Y_f);
B2j_R2  = jacobian(System_of_Eq, Y_b);
C1j_R2 = jacobian(System_of_Eq, Epsilon_t);

%=====================
common_params_SW;

%% regime 1 (PMP/AFP) 
disp('------------------------------')
disp( 'regime 1 (PMP/AFP) '   )
disp('------------------------------')
% set parameter values to matrix 

% p_st(:,1)= p_st(:,i);
 Regime_parameters_R1;

B1{1} = eval(B1j_R1);
A1{1} = eval(A1j_R1);
B2{1} = eval(B2j_R1);
C1{1} = eval(C1j_R1);

%                 A=B1^(-1)*A1
%                 B=B1^(-1)*B2
%                 C=B1^(-1)*C1
i=1;
            A{i,1}=B1{i}\A1{i};
            B{i,1}=B1{i}\B2{i};
            C{i,1}=B1{i}\C1{i};  

%% regime 2 (AMP/PFP) 
disp('------------------------------')
disp( 'regime 2 (AMP/PFP) '   )
disp('------------------------------')
% set parameter values to matrix 

%parameter_R1;
Regime_parameters_R2;
% Cap_phi_g1 = Cap_phi/(1+gam_st);            
 
 B1{2} = eval(B1j_R2);
A1{2} = eval(A1j_R2);
B2{2} = eval(B2j_R2);
C1{2} = eval(C1j_R2);

%                 A=B1^(-1)*A1
%                 B=B1^(-1)*B2
%                 C=B1^(-1)*C1
i=2;
            A{i,1}=B1{i}\A1{i};
            B{i,1}=B1{i}\B2{i};
            C{i,1}=B1{i}\C1{i};  

%% MOD method under MSRE using FM

  %  transtion Probs
  MS_transition_prob;
 
Opt.IRT=50;  %  Horizon of IRF
 
        [DET,FCC,OmegaK,GammaK,FK,~,IRF, IRF_PathWise]=fmmsre2(PP, A, B, C, [], Opt);  
    % (1) Area for DET                            
                    if DET(1,1)<1+1e-5 && DET(1,2)<=1-1e-5,                         
                         disp('     (1)  Solution is DET'),disp(DET)    
                         disp([ ' , since ', ...
                             num2str(DET(1,1)) '<1 and ' num2str(DET(1,2)) '<=1'     ]);   disp(' ');         end                    
                % (2) Area for INDET                    
                    if DET(1,1)<1+1e-5 && DET(1,2)>1-1e-5,                         
                         disp('     (2)  Solution is  INDET'),disp(DET)   ; 
                          disp(' ');          end
                % (3) Area for NSS                               
                    if DET(1,1) >1+1e-5 ,                         
                         disp('     (3)  Solution is   NSS'),disp(DET) ;   
                             disp(' ');        end        

%% The MOD Method under Fixed Regime 
      % regime 1
            [DET_R1,FCC_R1,OmegaK_R1,~,~,IRF_R1]=fmlre(A{1,1},B{1,1},C{1,1}, [], Opt);  
            if isnan(DET_R1(1,1)), [DET_R1,~,~,~,Geig_R1]=qzmlre(A{1,1},B{1,1},C{1,1}); end
      % regime 2      
            [DET_R2,FCC_R2,OmegaK_R2,~,~,IRF_R2]=fmlre(A{2,1},B{2,1},C{2,1}, [], Opt);  
            if isnan(DET_R2(1,1)), [DET_R2,~,~,~,Geig_R2]=qzmlre(A{2,1},B{2,1},C{2,1}); end

% (1) Area for DET       

                    if (DET_R1(1,1) < 1+1e-5 ) && (DET_R1(1,2)< 1-1e-5)                         
                         disp('     (1)  Solution is DET_R1')    ;  disp(DET_R1);
                                end                    
                % (2) Area for INDET_R1                    
                    if DET_R1(1,1)<1+1e-5 && DET_R1(1,2)>1-1e-5;                         
                         disp('     (2)  Solution is  INDET_R1'),disp(DET_R1)   ; 
                          disp(' ');          end
                % (3) Area for NSS                               
                    if DET_R1(1,1) >= 1+1e-5;                         
                         disp('     (3)  Solution is   NSS'),disp(DET_R1) ;   
                             disp(' ');        end

% (1) Area for DET                            
                    if DET_R2(1,1)<1+1e-5 && DET_R2(1,2)<=1-1e-5,                         
                         disp('     (1)  Solution is DET_R2'),disp(DET_R2)    
                         disp([ ' , since ', ...
                             num2str(DET_R2(1,1)) '<1 and ' num2str(DET_R2(1,2)) '<=1'     ]);   disp(' ');         end                    
                % (2) Area for INDET_R2                    
                    if DET_R2(1,1)<1+1e-5 && DET_R2(1,2)>1-1e-5,                         
                         disp('     (2)  Solution is  INDET_R2'),disp(DET_R2)   ; 
                          disp(' ');          end
                % (3) Area for NSS                               
                    if DET_R2(1,1) >1+1e-5 ,                         
                         disp('     (3)  Solution is   NSS'),disp(DET_R2) ;   
                             disp(' ');        end
                         
%% Plot  IRF
       %     mc          ${\mu_p}$       (long_name='gross price markup') 
%     zcap        ${z}$           (long_name='Capital utilization rate') 
%     rk          ${r^{k}}$       (long_name='rental rate of capital') 
%     k           ${k^{s}}$       (long_name='Capital services') 
%     pk          ${q}$           (long_name='real value of existing capital stock') 
%     c           ${c}$           (long_name='Consumption')
%     inve        ${i}$           (long_name='Investment')
%     y           ${y}$           (long_name='Output')
%     lab         ${l}$           (long_name='hours worked')
%     pinf        ${\pi}$         (long_name='Inflation')
%     w           ${w}$           (long_name='real wage')
%     r           ${r}$           (long_name='nominal interest rate')
%       Y_t = transpose([ mc zcap rk k pk c inve y lab pinf w r kp ]);  
       var_list ={  'gross price markup', 'Capital utilization rate', 'rental rate of capital','Capital services',...
                          'real value of existing capital stock', 'Consumption', 'Investment', 'Output',...
                          'hours worked','Inflation','real wage', 'nominal interest rate'} %, 'Capital stock'};
                      
  % varexo ea       ${\eta^a}$      (long_name='productivity shock')
%     eb          ${\eta^b}$      (long_name='preference shock')
%     eg          ${\eta^g}$      (long_name='Spending shock')
%     eqs         ${\eta^i}$      (long_name='Investment-specific technology shock')
%     em          ${\eta^m}$      (long_name='Monetary policy shock')
%     epinf       ${\eta^{p}}$    (long_name='Price markup shock')  
%     ew          ${\eta^{w}}$    (long_name='Wage markup shock')                      
 %      Epsilon_t = transpose([  a b g qs ms   spinf sw ]);
       shock_list ={'productivity shock' ,'preference shock', 'Spending shock',...
                           'Investment-specific technology shock','MP shock','Price markup shock','Wage markup shock'};
       nvar = length(var_list);
       nshock= length(shock_list);

    for j = 1:length(shock_list)
       figure('Name','IRF','FileName',[shock_list{j} ],'Position',[ 100, 100, 1200,600 ] )

        for i =1:length(var_list)
         subplot(4,3,i)
            plot(IRF{1,1}(:,(j-1)+1+(i-1)*nshock),'b-','LineWidth',2); hold on
            plot(IRF{2,1}(:,(j-1)+1+(i-1)*nshock),'r--','LineWidth',2); 
            xlim([0 Opt.IRT])

            if  i== nvar 
                 legend( { [ 'Regime AM: P= ' num2str(PP(1,1))   ], ...
                                [ 'Regime PM: P='  num2str(PP(2,2))   ] }  , 'Box', 'off',"FontSize",11)
            end
            title( [ var_list{ i } ] )
            set(gca,"FontSize",13);
            hold off
        end
              sgtitle(shock_list{j});
        exportgraphics(gcf,char(strcat('./graph/IRF_MS_',shock_list{j}, '.pdf')),'BackgroundColor','none')

    end    

%% Path-wise IRF    
    for j = 1:length(shock_list)
       figure('Name','Path-wise IRF','FileName',[shock_list{j} ],'Position',[ 100, 100, 1200,600 ] )

        for i =1:length(var_list)
         subplot(4,3,i)
            plot(IRF_PathWise{1,1}(:,(j-1)+1+(i-1)*nshock),'b-','LineWidth',2); hold on
            plot(IRF_PathWise{2,1}(:,(j-1)+1+(i-1)*nshock),'r--','LineWidth',2); 
            xlim([0 Opt.IRT])

            if  i== nvar
                 legend( { [ 'Regime AM: P= ' num2str(PP(1,1))   ], ...
                                [ 'Regime PM: P='  num2str(PP(2,2))   ] }  , 'Box', 'off',"FontSize",11)
            end
            title( [ var_list{ i } ] )
            set(gca,"FontSize",13);
            hold off
        end
              sgtitle(shock_list{j});
        exportgraphics(gcf,char(strcat('./graph/IRF_MS_PathWise_',shock_list{j}, '.pdf')),'BackgroundColor','none')

    end    


   %% Plot  IRF: Fixed Regime 
    for j = 1:length(shock_list)
       figure('Name','IRF','FileName',[shock_list{j} ],'Position',[ 100, 100, 1200,600 ] )

        for i =1:length(var_list)
         subplot(4,3,i)
            plot(IRF_R1(:,(j-1)+1+(i-1)*nshock),'b-','LineWidth',2); hold on
            plot(IRF_R2(:,(j-1)+1+(i-1)*nshock),'r--','LineWidth',2); 
            xlim([0 Opt.IRT])

            if  i== nvar 
                 legend( { [ 'Fix: Regime AM'   ], ...
                                [ 'Fix: Regime  PM' ] }  , 'Box', 'off',"FontSize",11)
            end
            title( [ var_list{ i } ] )
            set(gca,"FontSize",13);
            hold off
        end
              sgtitle(shock_list{j});
        exportgraphics(gcf,char(strcat('./graph/IRF_Fixed_',shock_list{j}, '.pdf')),'BackgroundColor','none')

    end    

                         
                         
                         