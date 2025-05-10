clear all
close all

%addpath('./Func/')
warning('off')

addpath(genpath('./MODMethod'))


   %%
   for  Panel=1:2 % choose from 1 to 6

         if Panel==1,        theta1=1.5; theta2=0.9;  theta3=0.0593; theta4=0.0593; 
                                     x_i=1; y_i=2;  
                                   disp('Panel A:  psi_pi (1) vs psi_pi (2) ')
         elseif Panel==2,     theta1=1.5; theta2=1.0;  theta3=0.0593; theta4=0.0593;
                                    x_i=1; y_i=3;                        
                                 disp('Panel B:  psi_pi (1)  vs delta b (1)  ')    

         end 

         %%
   
    %%
    vn{1}='\psi \pi (1)'; vn{2}='\psi \pi (2)'; vn{3}='\delta b (1)'; vn{4}='\delta b (2)';  % Parameter names
    vn{5}='\phi_B'; vn{6}='\phi_y';  vn{7}='\sigma'; 

        par=[theta1 theta2 theta3 theta4       ];  % Parameter bvlues
        par_L =[ 0.5     0.5    0.00     0.0     ];   % Lower bound
        par_U=[ 2.0     2.0    0.1       0.1       ];   % Upper bound
        disp('    Result with the following parameter values for example')
        disp(vn), disp(par)

%         
         [P, A,B,C]=func_Model_SW(par);
        [DET,FCC,OmegaK,GammaK]=X_SOLV_Q(P, A,B,C);      
       % (1) Area for DET                            
                    if DET(1,1)<=1 && DET(1,2)<=1,                         
                         disp('     (1)  Solution is DET'),disp(DET)    
                         disp([ ' , since ', ...
                             num2str(DET(1,1)) '<1 and ' num2str(DET(1,2)) '<=1'     ]);   disp(' ');         end                    
                % (2) Area for INDET                    
                    if DET(1,1)<=1 && DET(1,2)>1,                         
                         disp('     (2)  Solution is  INDET'),disp(DET)   ; 
                          disp(' ');          end
                % (3) Area for NSS                               
                    if DET(1,1) >1 ,                         
                         disp('     (3)  Solution is   NSS'),disp(DET) ;   
                             disp(' ');        end

        
        % Bounds for regime-switching parameters : 
        x_L=par_L(x_i); x_U=par_U(x_i);
        y_L=par_L(y_i); y_U=par_U(y_i);
        if  Panel==1
                xstep=0.05; ystep=0.05; 
        else
                xstep=0.05; ystep=0.025; 
        end    
        %xstep=0.05; ystep=0.05;        
        x_N=1;  y_N=1;
         
        xL=(x_L:xstep:x_U)';    xT=size(xL,1);      
        yL=(y_L:ystep:y_U)';    yT=size(yL,1);     

        LRTP1=zeros(xT,yT);     LRTP2=zeros(xT,yT); 
        DETA1=zeros(xT,yT);     DETA2=zeros(xT,yT); 
    
        %
        tx=1:xT;   ty=1:yT;
        [xGrid, yGrid]=ndgrid(tx,ty);
        Griddim=size(xGrid);

        DETA1_temp=zeros(Griddim);     DETA2_temp=zeros(Griddim); 
        DETA1_R1_temp=zeros(Griddim);     DETA2_R1_temp=zeros(Griddim); 
        DETA1_R2_temp=zeros(Griddim);     DETA2_R2_temp=zeros(Griddim);

        %
        disp('  Start calculation of Region of solutions')

        parfor i = 1: numel(xGrid)

            par_temp=par;
            % Run forward method xT*yT times
            x=xL(xGrid(i));   y=yL(yGrid(i));      
            par_temp(x_i)=x; par_temp(y_i)=y;

%          EvalModelQ;
              [P, A,B,C]=func_Model_SW(par_temp);
              [DET,FCC,OmegaK,GammaK,~,DET_R1,~,DET_R2]=X_SOLV_Q(P,A,B,C );
            
            %% For DET of MS model
            if ~isnan(DET(1,1))
              DETA1_temp(i)=DET(1,1);  DETA2_temp(i)=DET(1,2);
%              FCCA1(xGrid(i),yGrid(i))=FCC(1,1);  FCCA2(xGrid(i),yGrid(i))=FCC(1,2);
             else
                DETA1_temp(i)=10; DETA2_temp(i)=10
            end        

           %%  For DET of Fixed regime
            if ~isnan(DET_R1(1,1))
               DETA1_R1_temp(i)=(DET_R1(1,1)>1+1e-5);  DETA2_R1_temp(i)=(DET_R1(1,2)>1-1e-5);
%              FCCA1(xGrid(i),yGrid(i))=FCC(1,1);  FCCA2(xGrid(i),yGrid(i))=FCC(1,2);
             else
                DETA1_R1_temp(i)=10; DETA2_R1_temp(i)=10
            end     
            %%
            if ~isnan(DET_R2(1,1))
               DETA1_R2_temp(i)=(DET_R2(1,1)>1+1e-5);  DETA2_R2_temp(i)=(DET_R2(1,2)>1-1e-5);
%              FCCA1(xGrid(i),yGrid(i))=FCC(1,1);  FCCA2(xGrid(i),yGrid(i))=FCC(1,2);
             else
                DETA1_R2_temp(i)=10; DETA2_R2_temp(i)=10
            end     

        end
       
        %%
        DETA1= DETA1_temp; % NSS = DARK GRAY
        DETA2= DETA2_temp; % INDET=Light Gray
        %======================
        DETA1_Fix = DETA1_R1_temp +  DETA1_R2_temp  ;  % NSS = DARK GRAY_
        DETA2_Fix = DETA2_R1_temp  + DETA2_R2_temp; % INDET=Light Gray
     

    %% Generate Figure 1

    figure('Name','Region MS Model','FileName', ['Region-MS-panel' num2str(Panel) ] )

        opt.maxK=3000;  
        [X,Y] = meshgrid(xL,yL);
        %if Panel==1, [CC,hh]=contourf(X,Y,2*LRTP2',2*[1 1]); end
            clear map;
            map=[0.7 0.4 0.8 0.3]'*ones(1,3);   
         %             
            colormap(map) 
               contourf(X,Y,1*DETA2',1*[1-1e-5 1-1e-5],'LineWidth',2); hold on,... .% INDET=Light Gray
%              colormap(hot(8))
               contourf(X,Y,2*DETA1',2*[1+1e-5 1+1e-5],'LineWidth',2);  % NSS = DARK GRAY
             hold on              
        
            plot([x_L x_U]',[y_N y_N]','k','LineWidth',0.5); hold on,...
            plot([x_N x_N]',[y_L y_U]','k','LineWidth',0.5); hold on,...   
            yline(0,'k','LineWidth',0.5);
            if Panel==1,
                     title([ 'Taylor rule \psi_{\pi}  (1) vs.  Taylor rule \psi_{\pi} (2) ' ])
            elseif Panel ==2
                     title([ 'Taylor rule \psi_{\pi}  (1) vs.  Bohn rule \delta_b (1) ' ])
            
            else
                     title( [ 'x axis: ' vn{x_i} ' vs y axis '  vn{y_i}    ] );
            end    
            xlim([x_L x_U+0.0001]), ylim([y_L y_U]) 
            xlabel(vn{x_i}),     ylabel(vn{y_i})  
            pause(1)
            set(gca,'FontSize',14)
       
    hold off
    pause(0.3);
     
   saveas(gca, [ './graph/Region-MS-panel' num2str(Panel)  '.fig'  ] );
  
   %%

    figure('Name','Region Fix Regime','FileName', ['Region-Fix-panel' num2str(Panel) ] )

        opt.maxK=3000;  
        [X,Y] = meshgrid(xL,yL);
        %if Panel==1, [CC,hh]=contourf(X,Y,2*LRTP2',2*[1 1]); end
            clear map;
            map=[0.7 0.4 0.8 0.3]'*ones(1,3);   
         %             
            colormap(map) 
               contourf(X,Y,1*DETA2_Fix',1*[1-1e-5 1-1e-5],'LineWidth',2); hold on,... .% INDET=Light Gray
%              colormap(hot(8))
               contourf(X,Y,2*DETA1_Fix',2*[1+1e-5 1+1e-5],'LineWidth',2);  % NSS = DARK GRAY
             hold on              
        
            plot([x_L x_U]',[y_N y_N]','k','LineWidth',0.5); hold on,...
            plot([x_N x_N]',[y_L y_U]','k','LineWidth',0.5); hold on,...   
            yline(0,'k','LineWidth',0.5);
            if Panel==1,
                     title([ 'Taylor rule \psi_{\pi}  (1) vs.  Taylor rule \psi_{\pi} (2) ' ])
            elseif Panel ==2
                     title([ 'Taylor rule \psi_{\pi}  (1) vs.  Bohn rule \delta_b (1) ' ])
            
            else
                     title( [ 'x axis: ' vn{x_i} ' vs y axis '  vn{y_i}    ] );
            end    
            xlim([x_L x_U+0.0001]), ylim([y_L y_U]) 
            xlabel(vn{x_i}),     ylabel(vn{y_i})  
            pause(1)
            set(gca,'FontSize',14)
       
    hold off
    pause(0.3);
       
   saveas(gca, [ './graph/Region-fix-panel' num2str(Panel)  '.fig'  ] );


   end
  


%% Objective Function
function [varargout]=X_SOLV_Q(P,A, B,C) 
%     global S n P 
        warning('off')
       [DET,FCC,OmegaK,GammaK,FK]=fmmsre(P,A,B,C); % Modified Forward Method for MSRE
    varargout(1)={DET};
    varargout(2)={FCC};

    % Forward Method under Fixed Regime 
    if nargout>2 
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