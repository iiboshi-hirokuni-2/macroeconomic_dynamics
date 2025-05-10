%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Sample code Implementing the MOD method for MSRE models
%% Acronyms
%   DA / DIA : Determinacy-Adimissible / Determinacy-Inadmissible
%   DET/INDET/NSS : Determinate / Indeterminate / No Stable Solution
%   FM / GB : Forward Method / Groebner Basis
%   FS / MOD : Forward Solution / MOD solution
%% Three different ways.
% (1) modmsre_FM_GB.m : FM and GB 
%   : A complete and efficient approach
%       1) Try FM first, and GB if necessary.
%       2) Try GB Only if FM fails to detect Classification (DET/INDET/NSS) of MSRE models 
% (2) fmmsre.m: FM only 
%   : Not complete, but Practical in all economic models because
%     - Virtually all economic models are DA and DET is by far the most
%       important characteristic of research interest.
%     - Failure of FM occurs only for a subset of DIA models.
% (3) gbmsre.m: GB only
%   : A complete but inefficient approach because it works only
%     for a very small set of MSRE models.
%     Refer to gbmsre_Testing_computation_time.m.
%% To run this code, you must do the following.
%   1) Install Singular Program and follow the instruction in modmsre_FM_GB.m
%   2) Select an Approach (See below) 
%   3) Select a sample "Model" (See below) 
%   4) Run this code and follow the instruction in matlab command prompt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% This code tests the three different Approaches using 7 examples
%% Partitioning MSRE models by 3 Criteria
%% into 7 mutually disjoint and exhaustive subsets         
%            Criterion                    Cases                  Model # Below 
%   (1) Determinacy-Admissibility       : DA / DIA                1 / 2
%   (2) Number of Stable Solutions      : DET/ INDET / NSS        1 / 2 / 3
%   (3) Real-Valuedness of MOD solution : Real / Complex-valued   1 / 2
%   *** Note 1: DA model cannot have a complex-valued MOD solution.
%   *** Note 2: DIA model cannot be determinate.
%% This code tests the three different approachs using 7 Atheoretical Models 
%   --Model--  -----Criteria-----    -----Classification by FM First-----   --Classification By GB -------
%                                    FS Exist?  FS=MOD?  Class identified?    GB required? Class identified?                                                                                    
%   Model=111 : DA   DET     Real    	Yes      Yes        Yes                 No           Yes      
%   Model=121 : DA   INDET   Real    	Yes      Yes        Yes                 No           Yes  
%   Model=131 : DA   NSS     Real    	Yes      Yes        Yes                 No           Yes  
%   Model=221 : DIA  INDET   Real     	Yes      unknown    Yes                 No           Yes   
%   Model=222 : DIA  INDET   Complex  	No       unknown    No                  Yes          Yes 
%   Model=231 : DIA  NSS     Real       Yes      unknown    No                  Yes          Yes        
%   Model=232 : DIA  NSS     Complex    No       unknown    No                  Yes          Yes 

clear
%% Selecting an Approach and a Model
    Approach=2; % Set 1 using modmsre_FM_GB.m. 2 using fmmsre.m, 3 using gbmsre.m
    Model=231;  % Set model number defined above
    
%% DET-Admissible examples 
if Model==111 || Model==121 || Model==131  
	n=2; S=2;  
    P11=0.95; P22=0.9;
    P=[P11  1-P11;1-P22 P22];
    A{1,1}=[0.75 -0.05;0.1 0.67]; A{2,1}=[0.6 0.05;0.1 0.5];
    A{1,2}=[0.7 -0.1;0.3 0.6]; A{2,2}=[0.7 0.1;0.1 0.7];
    B{1,1}=[0.2 -0.3;-0.1 0.3]; B{2,1}=[0.3 -0.21;0.2 0.4];

    if Model==111, alpha=1;
    elseif Model==121, alpha=0.9; 
    elseif Model==131, alpha=1.5; 
    end
    for Si=1:S,  B{Si,1}=alpha*B{Si,1};
       	for Sj=1:S, A{Si,Sj}=(1/alpha)*A{Si,Sj}; end   
    end
end

%% DET-Inadmissible examples : Real-Valued MOD solution
if Model==221 || Model==231  
    n=1; S=2;
   	P11=0.3061; P12=1-P11; P22=0.9; P21=1-P22; 
   	P=[P11 P12;P21 P22];  
   	A{1,1}=0.2553; A{2,1}=5.0568; A{1,2}=A{1,1}; A{2,2}=A{2,1};
   	B{1,1}=1.5225; B{2,1}=-0.5870;
        
  	DET=fmmsre(P,A,B);
   	if Model==231, alpha=(1/DET(1));
      	for Si=1:S,  B{Si,1}=alpha*B{Si,1}; 
           	for Sj=1:S, A{Si,Sj}=(1/alpha)*A{Si,Sj}; end
        end
 	end 
end

%% DET-Inadmissible examples : Complex-Valued MOD solution
%
if Model==222 || Model==232
    n=2; S=2;  
    P11=0.95; P22=0.9;
    P=[P11  1-P11;1-P22 P22];
    A{1,1}=[0.3 0.2;0.1 0.45]; A{2,1}=[0.4 -0.1;0.1 0.35];  
    A{1,2}=A{1,1}; A{2,2}=A{2,1};
    B{1,1}=[0.7 0.1;0.05 0.65]; B{2,1}=[0.8 -0.05;0.2 0.75];
    
        if Model==222, alpha=0.5;
            for Si=1:S,  B{Si,1}=alpha*B{Si,1}; 
                for Sj=1:S, A{Si,Sj}=(1/alpha)*A{Si,Sj}; end
            end
        end 
end

%% Three Approaches
    if Approach==1
        [DETC,FCC,OmegaK,FK,DETCMOD,OmegaMOD,FMOD,DETC_All,AllOmegas]=modmsre_FM_GB(P,A,B); 
    elseif Approach==2
        [DETC,FCC,OmegaK,GammaK,FK,OtherOutput,IRS]=fmmsre(P,A,B);
        disp(' DETC information computed by Forward Method.')
        disp(DETC)
        
    elseif Approach==3
        [DETCMOD,OmegaMOD,FMOD,DET_All,AllOmegas]=gbmsre(P,A,B); 
    end

