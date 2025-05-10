clear all
close all

%addpath('./Func/')
warning('off')

addpath(genpath('./MODMethod'))

folder_name = "./graph";
if not(exist(folder_name,'dir'))
    mkdir(folder_name)
end


x_N=1;      
y_N=1;
xT = 20;   
yT = 20;
tx=1:xT;   
ty=1:yT;

for  Panel=1:3 % choose from 1 to 6
     if Panel==1       
         theta1=1.50; theta2=0.00;  theta3=0.05; theta4=0.00;  
         x_i=1; y_i=2;  
         disp('Panel A:  psi_pi (1) vs psi_pi (2) ')
     elseif Panel==2     
         theta1=1.50; theta2=0.00;  theta3=0.05; theta4=0.00;  
         x_i=3; y_i=4;           
         disp('Panel B:  delta_b (1)  vs delta_b (2)  ')
     elseif Panel==3     
         theta1=1.50; theta2=0.00;  theta3=0.05; theta4=0.00;  
         x_i=1; y_i=3;           
         disp('Panel C:  psi_pi (1)  vs delta b (1)  ')
     end 

    %
    vn{1}='\psi \pi (1)'; vn{2}='\psi \pi (2)'; vn{3}='\delta b (1)'; vn{4}='\delta b (2)';  % Parameter names


    par=[theta1 theta2 theta3 theta4];  % Parameter bvlues
    par_L=[0.00    0.00    0.00     0.00];   % Lower bound
    par_U=[2.25    1.05    0.10     0.11];   % Upper bound
    disp('    Result with the following parameter values for benchmark')
    disp(vn), disp(par)

      
    NK_model_with_habit;
    MS_transition_prob;
    common_params;

    B1j_R1 = -1*jacobian(System_of_Eq, Y_t);
    A1j_R1 = jacobian(System_of_Eq, Y_f);
    B2j_R1  = jacobian(System_of_Eq, Y_b);  
    C1j_R1 = jacobian(System_of_Eq, Epsilon_t);

    [A,B,C]=func_SetParam(par, B1j_R1,  A1j_R1,  B2j_R1,  C1j_R1);
    [DET,FCC,OmegaK,GammaK,FK,DetR1,FCCR1,DetR2,FCCR2]=X_SOLV_Q(PP, A,B,C); 
    % ----------   Benchmark parameters -------------------------------
    disp('     Fixed-Regimed   ')
    disp('     Regime 1: ')
    if DetR1(1,1)<1 && DetR1(1,2)<=1                         
         disp('     (1)  Solution is DET'),disp(DetR1)    
         disp([ ' , since ', ...
         num2str(DetR1(1,1)) '<1 and ' num2str(DetR1(1,2)) '<=1'     ]);   disp(' ');        
    end                    
    % (2) Area for INDET                    
    if DetR1(1,1)<1 && DetR1(1,2)>1                         
         disp('     (2)  Solution is  INDET'),disp(DetR1)   ; 
         disp(' ');          
    end
    % (3) Area for NSS                               
    if DetR1(1,1) >1                          
         disp('     (3)  Solution is   NSS'),disp(DetR1) ;       
         disp(' ');       
    end

    disp('     Regime 2: ')
    if DetR2(1,1)<1 && DetR2(1,2)<=1                         
         disp('     (1)  Solution is DET'),disp(DetR2)    
         disp([ ' , since ', ...
         num2str(DetR2(1,1)) '<1 and ' num2str(DetR2(1,2)) '<=1'     ]);   disp(' ');        
    end                    
    % (2) Area for INDET                    
    if DetR2(1,1)<1 && DetR2(1,2)>1                         
         disp('     (2)  Solution is  INDET'),disp(DetR2)   ; 
         disp(' ');          
    end
    % (3) Area for NSS                               
    if DetR2(1,1) >1                          
         disp('     (3)  Solution is   NSS'),disp(DetR2) ;       
         disp(' ');       
    end



    disp('     Markov-Switching   ')
    % (1) Area for DET                            
    if DET(1,1)<1 && DET(1,2)<=1                         
         disp('     (1)  Solution is DET'),disp(DET)    
         disp([ ' , since ', ...
         num2str(DET(1,1)) '<1 and ' num2str(DET(1,2)) '<=1'     ]);   disp(' ');        
    end                    
    % (2) Area for INDET                    
    if DET(1,1)<1 && DET(1,2)>1                         
         disp('     (2)  Solution is  INDET'),disp(DET)   ; 
         disp(' ');          
    end
    % (3) Area for NSS                               
    if DET(1,1) >1                          
         disp('     (3)  Solution is   NSS'),disp(DET) ;       
         disp(' ');       
    end



    % Bounds for regime-switching parameters : 
    x_L=par_L(x_i); x_U=par_U(x_i);
    y_L=par_L(y_i); y_U=par_U(y_i);

    xL = transpose(linspace(x_L,x_U,xT));
    yL = transpose(linspace(y_L,y_U,yT));

    [xGrid, yGrid]=ndgrid(tx,ty);
    Griddim=size(xGrid);

    DETA1_temp=zeros(Griddim);     
    DETA2_temp=zeros(Griddim); 
    DETA1_temp_R1=zeros(Griddim);     
    DETA2_temp_R1=zeros(Griddim); 

    disp('  Start calculation of Region of solutions')

    parfor i = 1: numel(xGrid)
        par_temp=par;
        % Run forward method xT*yT times
        x=xL(xGrid(i));   
        y=yL(yGrid(i));      
        par_temp(x_i)=x; 
        par_temp(y_i)=y;

%       EvalModelQ;
        [A,B,C]=func_SetParam(par_temp, B1j_R1,  A1j_R1,  B2j_R1,  C1j_R1);
        [DET,FCC,OmegaK,GammaK,FK,DetR1,FCCR1,DetR2,FCCR2]=X_SOLV_Q(PP,A,B,C );

        % For DET 
        if ~isnan(DET(1,1))
            DETA1_temp(i)=DET(1,1);  DETA2_temp(i)=DET(1,2);
%              FCCA1(xGrid(i),yGrid(i))=FCC(1,1);  FCCA2(xGrid(i),yGrid(i))=FCC(1,2);
        else
            DETA1_temp(i)=10; DETA2_temp(i)=10;
        end

        if Panel == 3
            if ~isnan(DetR1(1,1))
                DETA1_temp_R1(i)=DetR1(1,1);  DETA2_temp_R1(i)=DetR2(1,2);
%              FCCA1(xGrid(i),yGrid(i))=FCC(1,1);  FCCA2(xGrid(i),yGrid(i))=FCC(1,2);
            else
                DETA1_temp_R1(i)=10; DETA2_temp_R1(i)=10;
            end
        end
                     
    end
    DETA1= DETA1_temp;   
    DETA2= DETA2_temp; 
    DETA1_R1= DETA1_temp_R1;   
    DETA2_R1= DETA2_temp_R1; 
    No_Solution = 3.1*(DETA2_R1>=10);

    figure('Name','Region','FileName', ['Region-panel' num2str(Panel) ] )
    opt.maxK=3000;  
    [X,Y] = meshgrid(xL,yL);
    clear map;
    map=[0.7 0.4 0.8 0.3]'*ones(1,3); 
     
    colormap(map) 
    contourf(X,Y,1*DETA2',1*[1 1],'LineWidth',2); hold on,... .% INDET=Light Gray
    contourf(X,Y,2*DETA1',2*[1 1],'LineWidth',2);  % NSS = DARK GRAY
    hold on    
    
    plot([x_L x_U]',[y_N y_N]','k','LineWidth',0.5); hold on,...
    plot([x_N x_N]',[y_L y_U]','k','LineWidth',0.5); hold on,...   
    yline(0,'k','LineWidth',0.5);
    
    if Panel==1
        title([ 'Taylor rule \psi_{\pi}  (1) vs.  Taylor rule \psi_{\pi} (2) ' ]);
    elseif Panel ==2
        title([ 'Bohn rule \delta_b (1) vs.  Bohn rule \delta_b (2) ' ]);
    elseif Panel ==3    
        title([ 'Taylor rule \psi_{\pi}  (1) vs.  Bohn rule \delta_b (1) ' ]);
    else
        title( [ 'x axis: ' vn{x_i} ' vs y axis '  vn{y_i}    ] );
    end
    xlim([x_L x_U]), ylim([y_L y_U]) 
    xlabel(vn{x_i}),     ylabel(vn{y_i})
    %pause(1)
    set(gca,'FontSize',12)
    
    hold off
    %pause(0.3);    
    saveas(gca, [ './graph/Region-panel' num2str(Panel)  '.fig'  ] );

    if Panel == 3
        figure('Name','Region','FileName', ['Region-panel' num2str(Panel),'(Fixed)' ] )
        opt.maxK=3000;  
        [X,Y] = meshgrid(xL,yL);
        clear map;
        map=[0.7 0.4 0.8 0.3]'*ones(1,3); 
         
        colormap(map) 
        contourf(X,Y,1*DETA2_R1',1*[1 1],'LineWidth',2); hold on,... .% INDET=Light Gray
        contourf(X,Y,2*DETA1_R1',2*[1 1],'LineWidth',2);  % NSS = DARK GRAY
        hold on    
        %contourf(X,Y,No_Solution',2.4*[1 1],'LineWidth',2);  % No Solution = BLACK
        
        plot([x_L x_U]',[y_N y_N]','k','LineWidth',0.5); hold on,...
        plot([x_N x_N]',[y_L y_U]','k','LineWidth',0.5); hold on,...   
        yline(0,'k','LineWidth',0.5);

        title([ 'Taylor rule \psi_{\pi}  (1) vs.  Bohn rule \delta_b (1) ' ]);

        xlim([x_L x_U]); ylim([y_L y_U]);
        xlabel(vn{x_i});     ylabel(vn{y_i});
        %pause(1)
        set(gca,'FontSize',12)
    
        hold off
        %pause(0.3);  
        saveas(gca, [ './graph/Region-panel_Fixed' num2str(Panel)  '.fig'  ] );
        
    end
end





%% Objective Function
function [varargout]=X_SOLV_Q(P,A, B,C) 
%     global S n P 
        warning('off')
       [DET,FCC,OmegaK,GammaK,FK]=fmmsre(P,A,B,C); % Modified Forward Method for MSRE
    varargout(1)={DET};
    varargout(2)={FCC};

    % Forward Method under Fixed Regime 
    if nargout>5 
        % Regime 1
            A_R1=A{1,1}; B_R1=B{1,1}; C_R1=C{1,1}; 
            [DET_R1,FCC_R1]=fmlre(A_R1,B_R1,C_R1);  % Modified Forward Method for LRE
        % Regime 2  
            A_R2=A{2,1}; B_R2=B{2,1}; C_R2=C{1,1};
            [DET_R2,FCC_R2]=fmlre(A_R2,B_R2,C_R2);  % Modified Forward Method for LRE
         
    % Output         
        varargout(3)={OmegaK};
        varargout(4)={GammaK};
        varargout(5)={FK};      
        varargout(6)={DET_R1};
        varargout(7)={FCC_R1};
        varargout(8)={DET_R2}; 
        varargout(9)={FCC_R2};  
        varargout(10)={A};  
        varargout(11)={B};  
    end
end