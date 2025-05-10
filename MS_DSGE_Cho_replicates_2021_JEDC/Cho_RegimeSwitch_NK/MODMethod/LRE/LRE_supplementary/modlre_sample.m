%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Sample code Implementing the MOD method for LRE models
%% Acronyms
%   DA / DIA : Determinacy-Adimissible / Determinacy-Inadmissible
%   DET/INDET/NSS : Determinate / Indeterminate / No Stable Solution
%   FM / QZ : Forward Method / QZ Method
%   FS / MOD : Forward Solution / MOD solution
%% Three different Approaches to MOD method for LRE models.
% (1) modlre_FM_QZ.m : FM and QZ 
%   : A complete and efficient approach
%       1) Try FM first, and QZ if necessary.
%       2) Try QZ Only if FM fails to detect Classification (DET/INDET/NSS) of LRE models 
% (2) fmlre.m: FM only 
%   : Not complete, but as efficient as QZ method in all economic models because
%     - Virtually all economic models are DA and DET is by far the most
%       important characteristic of research interest.
%     - Failure of FM occurs only for a subset of DIA models.
% (3) qzmlre.m: QZ only
%   : A complete and efficient stand-alone approach 
%     comparable to gensys but fits in our model and solution representation
%% To run this code, you must do the following.
%   2) Select an Approach (See below) 
%   3) Select a sample "Model" (See below) 
%   4) Run this code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% This code tests the three different Approaches using 8 examples
%   ** Two DA-DET examples are tested to use the same models in modlre_methods_comparison.m 
%% Partitioning LRE models by 3 Criteria
%% into 7 mutually disjoint and exhaustive subsets         
%            Criterion                    Cases                  Model # Below 
%   (1) Determinacy-Admissibility       : DA / DIA                1 / 2
%   (2) Number of Stable Solutions      : DET/ INDET / NSS        1 / 2 / 3
%   (3) Real-Valuedness of MOD solution : Real / Complex-valued   1 / 2
%   *** Note 1: DA model cannot have a complex-valued MOD solution.
%   *** Note 2: DIA model cannot be determinate.
%% This code tests the three different approachs using 7 Atheoretical Models 
%   --Model--  -----Criteria-----    -----Classification by FM First-----   --Classification By QZ -----       --modlre_methods_comparison.m--
%                                    FS Exist?  FS=MOD?  Class identified?    QZ required? Class identified?                                                                                    
%   Model=1111: DA   DET     Real    	Yes      Yes        Yes                 No           Yes                 Case 1        
%   Model=1112: DA   DET     Real    	Yes      Yes        Yes                 No           Yes                 Case 2       
%   Model=121 : DA   INDET   Real    	Yes      Yes        Yes                 No           Yes                 Case 1 transformed to be indeterminate
%   Model=131 : DA   NSS     Real    	Yes      Yes        Yes                 No           Yes                 Case 1 transformed to have NSS
%   Model=221 : DIA  INDET   Real     	Yes      unknown    Yes                 No           Yes                 Case 3 transformed to be indeterminate  
%   Model=222 : DIA  INDET   Complex  	No       unknown    No                  Yes          Yes                 Case 4 transformed to have NSS
%   Model=231 : DIA  NSS     Real       Yes      unknown    No                  Yes          Yes                 Case 3
%   Model=232 : DIA  NSS     Complex    No       unknown    No                  Yes          Yes                 Case 4

clear
%% Selecting an Approach and a Model
    Approach=1; % Set 1 using modlre_FM_QZ.m. 2 using fmmsre.m, 3 using qzmlre.m
    Model=1111;  % Set model number defined above

%% Models
%  DET-Admissible examples 
    if Model==1111 || Model==1112 || Model==121 || Model==131  
        if Model==1111, B0=[1 0       ;0.1   0.9];       A0=[0 0;0 1]; B1=[1.1 0.5;0 0];  alpha=1; end
        if Model==1112, B0=[1 0.000001;0   0.9];        A0=[0 0;0 1]; B1=[1.1 0;0 0]; alpha=1;  end
        if Model==121, B0=[1 0       ;0.1   0.9];       A0=[0 0;0 1]; B1=[1.1 0.5;0 0];  alpha=1.5; end           
        if Model==131, B0=[1 0       ;0.1   0.9];       A0=[0 0;0 1]; B1=[1.1 0.5;0 0];  alpha=0.5; end    
    end
%  DET-Inadmissible examples : Real-Valued MOD solution
    if Model==221, B0=[1 0       ;0   0.9];         A0=[0 0;0 1]; B1=[1.1 0;0 0]; alpha=1.5;  end
    if Model==231, B0=[1 0       ;0   0.9];         A0=[0 0;0 1]; B1=[1.1 0;0 0]; alpha=1;  end
% DET-Inadmissible examples : Complex-Valued MOD solution
    if Model==222, a1=0.1; a2=-a1;
              B0=[1 0       ;a2    1];         A0=[0 0;0 1]; B1=[1 a1;0 0]; alpha=2;  end

    if Model==232, a1=0.1; a2=-a1;
              B0=[1 0       ;a2    1];         A0=[0 0;0 1]; B1=[1 a1;0 0]; alpha=1;  end
          
    A=(1/alpha)*B0\A0; B=alpha*B0\B1;       
%% Three Approaches
    if Approach==1
        [DETC,FCC,OmegaK,FK,DETCMOD,OmegaMOD,FMOD]=modlre_FM_QZ(A,B);
    elseif Approach==2
        [DETC,FCC,OmegaK,GammaK,FK]=fmlre(A,B);
        disp(' DETC information computed by Forward Method.')
        disp(DETC)
    elseif Approach==3
        [DETCMOD,OmegaMOD,GammaMOD,FMOD,Geig,gvAll,OmegaAll,GammaAll,FAll]=qzmlre(A,B); 
        disp(' DETCMOD information computed by QZ Method.')
        disp(DETCMOD)
    end
