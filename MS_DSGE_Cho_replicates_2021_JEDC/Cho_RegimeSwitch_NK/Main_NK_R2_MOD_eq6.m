clc

% addpath('./TZMatlabPrograms/')
addpath(genpath('./MODMethod'))
addpath('./Func/')

clear all
 close all

disp( '  '  )
disp('============================================')
disp('           start  NK-DSGE model ')
disp('============================================')

%NK_model;
NK_model_with_habit;
%NK_model_with_habit_blanced_budget;
common_params;


% regime 1 (AMP/PFP) 
disp('------------------------------')
disp( 'regime 1 (AMP/PFP) '   )
disp('------------------------------')
% set parameter values to matrix 
B1j_R1 = -1*jacobian(System_of_Eq, Y_t);
A1j_R1 = jacobian(System_of_Eq, Y_f);
B2j_R1  = jacobian(System_of_Eq, Y_b);
C1j_R1 = jacobian(System_of_Eq, Epsilon_t);
Regime_parameters_R1;
% parameter_R2;
%Cap_phi_g1 = Cap_phi/(1+gam_st);


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

% regime 2 (PMP/AFP) 
disp('------------------------------')
disp( 'regime 2 (PMP/AFP) '   )
disp('------------------------------')
% set parameter values to matrix 

%parameter_R1;
B1j_R2 = -1*jacobian(System_of_Eq, Y_t);
A1j_R2 = jacobian(System_of_Eq, Y_f);
B2j_R2  = jacobian(System_of_Eq, Y_b);
C1j_R2 = jacobian(System_of_Eq, Epsilon_t);
 Regime_parameters_R2;
%Cap_phi_g1 = Cap_phi/(1+gam_st);            
 
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
 
Opt.IRT=80;  %  Horizon of IRF

        %[DET,FCC,OmegaK,GammaK,FK,~,IRF]=fmmsre(PP, A, B, C, [], Opt);  
        [DET,FCC,OmegaK,GammaK,FK,~,IRF, IRF_PathWise]=fmmsre2(PP, A, B, C, [], Opt);
init_state = 1;
T = 2500;
vecsig = p_er;
        [State_path,Endo_path,Exo_path] = stock_simul_RS(OmegaK,GammaK,PP,init_state,T,vecsig);

            

%% The MOD Method under Fixed Regime 
      % regime 1 
            [DET_R1,FCC_R1,OmegaK_R1,~,~,IRF_R1]=fmlre(A{1,1},B{1,1},C{1,1},[],Opt);  
            if isnan(DET_R1(1,1)), [DET_R1,~,~,~,Geig_R1]=qzmlre(A{1,1},A{1,1}); end
      % regime 2      
            [DET_R2,FCC_R2,OmegaK_R2,~,~,IRF_R2]=fmlre(A{2,1},B{2,1},C{2,1},[],Opt);  
            if isnan(DET_R2(1,1)), [DET_R2,~,~,~,Geig_R2]=qzmlre(A{2,1},B{2,1}); end

% [2] Indicators for DET/INDET/NSS and DET-Inadmissibility
%% : Decision making about Classification of MSRE models.
%   [1] Forward Method using this code.
%       (1) If DETC1*DETC2<1, the model is DA. FS=MOD. Done.
%       (2) If DETC1*DETC2>=1 and DETC1<1, model is indeterminate. FS=MOD is not confirmed.
%       (3) If DETC1*DETC2>=1 and DETC1>=1, FM cannot tell DET/INDET/NSS.
%       (4) If DETC1=NaN, the forward fail to converge. FM cannot tell DET/INDET/NSS.
%   [2] Apply QZ method using qzmlre.m in the case of (3) from [1]
%       : Identify the MOD solution by solving all MSV solutions, completing the MOD method. 

                % (1) Area for DET
                    if DET_R1(1,1)<1 && DET_R1(1,2)<=1,  
                        disp('     (1-1) DET under Fixed regime 1'),disp(DET_R1);   
                        disp( [' , since ', ...
                             num2str(DET_R1(1,1)) '<1 and ' num2str(DET_R1(1,2)) '<1'     ]);    disp(' ');  end
                    if DET_R2(1,1)<1 && DET_R2(1,2)<=1,   
                        disp('     (1-2) DET under Fixed regime 2'),disp(DET_R2);   
                        disp( [' , since ', ...
                             num2str(DET_R2(1,1)) '<1 and ' num2str(DET_R2(1,2)) '<1'     ]); disp(' ');  end              
                    if DET(1,1)<1 && DET(1,2)<=1,                         
                         disp('     (2)   DET under Markov-switching'),disp(DET)    
                         disp([ ' , since ', ...
                             num2str(DET(1,1)) '<1 and ' num2str(DET(1,2)) '<1'     ]);   disp(' ');         end
                    
                % (2) Area for INDET
                      if DET_R1(1,1)<1 && DET_R1(1,2)>1,  
                       disp('     (1-1) INDET under Fixed regime 1'),disp(DET_R1);   
                        disp(' '); end
                    if DET_R2(1,1)<1 && DET_R2(1,2)>1,   
                        disp('     (1-2) INDET under Fixed regime 2'),disp(DET_R2);   
                         disp(' ');  end              
                    if DET(1,1)<1 && DET(1,2)>1,                         
                         disp('     (2)   INDET under Markov-switching'),disp(DET)   ; 
                          disp(' ');          end
                % (3) Area for NSS
                    if DET_R1(1,1)>1 , 
                       disp('     (1-1) NSS under Fixed regime 1'),disp(DET_R1);   
                        disp(' ');  end
                    if DET_R2(1,1)>1 ,   
                        disp('     (1-2) NSS under Fixed regime 2'),disp(DET_R2);   
                          disp(' '); end              
                    if DET(1,1) >1 ,                         
                         disp('     (2)   NSS under Markov-switching'),disp(DET) ;   
                             disp(' ');        end


       
       % Y_t =([ y pai R b tau e  g a mupR tp ]);  
       var_list ={  'output', 'Output gap', 'Natural Y', 'Natural rate', 'pi', 'marginal cost', 'R', 'Debt-to-GDP', 'Tax', 'Consumption', 'dyObs', 'dcObs', 'piObs', 'rnObs', 'MPshock' ,'GOVshock', 'TFPshock', 'COSTshock'};
       %Epsilon_t = ([  eps_MP  eps_g eps_a  eps_tau  eps_e    eps_mupR eps_tp  ]);
       shock_list ={'MPshock' ,'GOVshock', 'TFPshock', 'COSTshock'};
       obs_list ={'dyObs', 'dcObs', 'piObs', 'rnObs'};
       nvar = length(var_list);
       nshock= length(shock_list);
       nobs = length(obs_list);

       Obs_path = zeros(nobs,T);
       for i =1:nobs
           Obs_path(i,:) = Endo_path(nvar - nshock+i,:);
       end

           Old = cd('.\graph'); 
        figure
        plot(State_path);
        ylim([0.9 2.1])
        title('Policy regime')
        exportgraphics(gcf,char(strcat('Stock_simul_Regimes.pdf')),'BackgroundColor','none')

        figure
        for i=1:nshock
            subplot(nshock,1,i)
            plot(Exo_path(i,:));
            title( [ shock_list{ i } ] )
        end
        exportgraphics(gcf,char(strcat('Stock_simul_shocks.pdf')),'BackgroundColor','none')

        figure
        temp = nvar - nshock- nobs;
        for i=1:nshock
            subplot(nshock,1,i)
            plot(Endo_path(temp+i,:));
            title( [ shock_list{ i } ] )
        end
        exportgraphics(gcf,char(strcat('Stock_simul_exo_variables.pdf')),'BackgroundColor','none')

        figure
        temp = nvar - nshock;
        for i=1:nobs
            subplot(nobs,1,i)
            plot(Endo_path(temp+i,:));
            title( [ obs_list{ i } ] )
        end
        exportgraphics(gcf,char(strcat('Stock_simul_observation.pdf')),'BackgroundColor','none')

        % Plot  IRF
 
    for j = 1:nshock
        figure('Name','IRF','FileName',[shock_list{j} ],'Position',[ 100, 100, 1200,600 ] )

        for i =1:nvar-nshock-nobs
            subplot(4,3,i)
            plot(IRF{1,1}(:,(j-1)+1+(i-1)*nshock),'b-','LineWidth',2); hold on
            plot(IRF{2,1}(:,(j-1)+1+(i-1)*nshock),'r--','LineWidth',2); 
            xlim([0 Opt.IRT])
            title( [ var_list{ i } ] )
            set(gca,"FontSize",13);
            hold off
        end
        subplot(4,3,nvar-nshock-nobs+1)
        i = nvar-nshock + j;
        plot(IRF{1,1}(:,(j-1)+1+(i-1)*nshock),'b-','LineWidth',2); hold on
        plot(IRF{2,1}(:,(j-1)+1+(i-1)*nshock),'r--','LineWidth',2); 

        legend( { [ 'Regime AM/PF: P= ' num2str(PP(1,1))   ], ...
        [ 'Regime PM/AF: P='  num2str(PP(2,2))   ] }  , 'Box', 'off',"FontSize",11)
        title( [ shock_list{ j } ] )
         set(gca,"FontSize",13);
        hold off
        sgtitle(shock_list{j});
        exportgraphics(gcf,char(strcat('IRF_MS_',shock_list{j}, '.pdf')),'BackgroundColor','none')
        saveas(gcf,strcat('IRF_MS_',shock_list{j}, '.png'))

    end    

%% Path-wise IRF    
    for j = 1:length(shock_list)
       figure('Name','Path-wise IRF','FileName',[shock_list{j} ],'Position',[ 100, 100, 1200,600 ] )

        for i =1:nvar-nshock-nobs
            subplot(4,3,i)
            plot(IRF_PathWise{1,1}(:,(j-1)+1+(i-1)*nshock),'b-','LineWidth',2); hold on
            plot(IRF_PathWise{2,1}(:,(j-1)+1+(i-1)*nshock),'r--','LineWidth',2); 
            xlim([0 Opt.IRT])
            title( [ var_list{ i } ] )
            set(gca,"FontSize",13);
            hold off
        end
        subplot(4,3,nvar-nshock-nobs+1)
        i = nvar-nshock + j;
        plot(IRF_PathWise{1,1}(:,(j-1)+1+(i-1)*nshock),'b-','LineWidth',2); hold on
        plot(IRF_PathWise{2,1}(:,(j-1)+1+(i-1)*nshock),'r--','LineWidth',2); 
        legend( { [ 'Regime AM/PF: P= ' num2str(PP(1,1))   ], ...
        [ 'Regime PM/AF: P='  num2str(PP(2,2))   ] }  , 'Box', 'off',"FontSize",11)
        title( [ shock_list{ j } ] )
         set(gca,"FontSize",13);
        hold off
        sgtitle(shock_list{j});
        exportgraphics(gcf,char(strcat('IRF_MS_PathWise_',shock_list{j}, '.pdf')),'BackgroundColor','none')
        saveas(gcf,strcat('IRF_MS_PathWise_',shock_list{j}, '.png'))

    end    


   % Plot  IRF: Fixed Regime 
    for j = 1:length(shock_list)
        figure('Name','IRF','FileName',[shock_list{j} ],'Position',[ 100, 100, 1200,600 ] )

        for i =1:nvar-nshock-nobs
            subplot(4,3,i)
            plot(IRF_R1(:,(j-1)+1+(i-1)*nshock),'b-','LineWidth',2); hold on
            plot(IRF_R2(:,(j-1)+1+(i-1)*nshock),'r--','LineWidth',2); 
            xlim([0 Opt.IRT])
            title( [ var_list{ i } ] )
            set(gca,"FontSize",13);
            hold off
        end
        subplot(4,3,nvar-nshock-nobs+1)
        i = nvar-nshock + j;
        plot(IRF_R1(:,(j-1)+1+(i-1)*nshock),'b-','LineWidth',2); hold on
        plot(IRF_R2(:,(j-1)+1+(i-1)*nshock),'r--','LineWidth',2); 
        legend( { [ 'Regime AM/PF (Fixed)'   ], ...
        ['Regime PM/AF (Fixed)'  ] }  , 'Box', 'off',"FontSize",11)
        title( [ shock_list{ j } ] )
         set(gca,"FontSize",13);
        hold off
        sgtitle(shock_list{j});
        exportgraphics(gcf,char(strcat('IRF_Fixed_',shock_list{j}, '.pdf')),'BackgroundColor','none')
        saveas(gcf,strcat('IRF_Fixed_',shock_list{j}, '.png'))

    end    
    cd(Old);


