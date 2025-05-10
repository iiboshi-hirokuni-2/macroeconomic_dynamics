clear, clc
%% Just Run this code.
disp('========================================================')
disp('The main purpose of Cho(2019) is a complete classification of MSRE')
disp('models into DET/INDET/NSS. The main contribution is to derive')
disp('necessary and sufficient conditions for each case using the information')
disp('of the most stable solution (MOD solution) only.')
disp('The proposed MOD method also provides an identification condition for')
disp('the MOD solutution and a straightforward technique to compute the MOD')
disp('solution without relying on the eigensystem, which is unknown.')
disp('The methodology must nest the LRE models as special case and the ')
disp('proposed method must be identical to the existing technique relying')
disp('on eigensystem.')
disp('This code shows the equivalence numerically.')
disp(' To do so, (1) we first introduce definitions and acronyms')
disp(' followed by (2) partitioning of the LRE models based on MOD solution')
disp(' for identification of the MOD solution.')
disp(' and (3) compare the solution techniques.')
disp('========================================================')
disp('[1] Acronyms / Definition')
disp('  (1) Model')
disp('      LRE/MSRE : Linear / Markov-switching Rational Expectations') 
disp('      (P,A,B) : the matrices defining models. Refer to Paper.')
disp('  (2) Solution')
disp('      (Omega,F) : matrices defining the full set of solutions')
disp('             in terms of Omega: MSV solutions) ')
disp('             its uniquely associated F: sunspot components.')
disp('          MOD (Mininum of modulus) solution: the most stable solution.')
disp('          Omegaj,Fj = j-th stable solution and its associated F')
disp('          OmegaMOD, FMOD = MOD solution and its associated F.')
disp('          Omega_ns: Omega with the n smallest Generalized Eigenvalues.')
disp('  (3) Model Partitioning criteria')
disp('      1) DA / DIA : Determinacy-Adimissible / Determinacy-Inadmissible.')
disp('      2) DET/INDET/NSS : Determinate / Indeterminate / No Stable Solution')
disp('         --------------------------------------------------------------')
disp('                        DA                    DIA       ')
disp('                 -----------------       ---------------')
disp('          DET    :    DA/DET        :     Impossible    ')
disp('          INDET  :    DA/INDET      :     DIA/INDET     ')
disp('          NSS    :    DA/NSS        :     DIA/INDET     ')
disp('         --------------------------------------------------------------')
disp('      *** Note 1: DA model cannot have a complex-valued MOD solution.')
disp('      *** Note 2: DIA model cannot be determinate.')
disp('                  A complex-valued MOD solution automatically implies DIA.')
disp('          So, there are in total 7 mutially disjoint and exhausitive subsets of LRE models.')    
disp('  (4) MOD information')
disp('      r( ) : Spectral Radius')
disp('      Geig : Generalized eigenvalues')
disp('      r1(MOD) = r(OmegaMOD), r2(MOD)=r(FMOD)')
disp('      r1(j)=r(Omegaj), r2(j)=r(Fj)')
disp('  (5) Solution Techniques ')
disp('      FM / QZ /gensys: Forward Method / QZ Method /gensys alrogithm')
disp('========================================================')
disp('The following two partitioning are equivalent to each other and ')
disp('Classification of these two into DET/INDET/NSS is equivalent to that of gensys')
disp('Refer to Paper for Proof.')
disp('[1] Partitioning LRE models by MOD solution and the Identifying Conditions')
disp('      --------------------------------------------------------------------------------------------')
disp('                       DA             :       DIA              :  Pulling Both Together           ')
disp('               =======================:========================:==================================')
disp('                r1(MOD)*r2(MOD)<1     : r1(MOD)*r2(MOD)>=1     :                                  ')
disp('                & Real MOD            : Or Complex OmegaMOD    :                                  ') 
disp('               -----------------------:------------------------:----------------------------------')
disp('       DET     r1(MOD)<1 & r2(MOD)<=1 :    Impossible          : Real MOD & r1(MOD)<1 & r2(MOD)<=1')
disp('       INDET               r2(MOD)>1  :    r1(MOD)<1           : r1(MOD)<1 and r2(MOD)>1          ')
disp('       NSS     r1(MOD)>=1             :    r1(MOD)>=1          : r1(MOD)>=1                       ')              
disp('      --------------------------------------------------------------------------------------------')
disp('  ')
disp('[2] Partitioning LRE models by Eigensystem')
disp('    Let (xi_1,...,xi_2n) be generalized eigenvalues.')
disp('    Omega_ns=Omega(xi_1,...,xi_n)')
disp('    If Omega_ns exists, OmegaMOD  = Omega_ns and r1(MOD)=|xi_n|.')
disp('    If not,             OmegaMOD ~= Omega_ns and r1(MOD)=|xi_{n+j}|>|xi_n| for some j=1,..,n. ') 
disp('    Let xi_(n+J) be the eigenvalue of Omega_MOD. Then J=0 if Omega_ns exists and J>0 otherwise.')
disp('      ----------------------------------------------------------------------------------------')
disp('                       DA             :        DIA             :  Pulling Both Together       ')
disp('               =======================:========================:==============================')
disp('               |xi_n|<|xi_{n+1}| &    : |xi_n| = |xi_{n+1}| Or :                              ')
disp('               J=0 (Omega_ns exist)   : J>0 (Omega_ns ~exist)  :                              ')
disp('               -----------------------:------------------------:------------------------------')
disp('       DET     |xi_n|<1<=|xi_{n+1}|   :   Impossible           : J=0 & |xi_n|<1<=|xi_{n+1}|   ')
disp('       INDET             |xi_{n+1}|<1 :   |xi_{n+J}|<1         : |xi_{n+J}|<1 and |xi_{n+1}|<1')
disp('       NSS     |xi_n|>=1              :   |xi_{n+J}|>=1        : |xi_{n+J}|>=1                ')
disp('      ----------------------------------------------------------------------------------------')
disp('  ')  
disp('[3] Classification by gensys')
disp('          DET if eu=[1 1],   INDET if eu=[1 0],  NSS if eu(1)<1.')
disp('  ')  
disp('     Note: 1) Existence of Omega_ns (J=0) is an additional condition for DET.')
disp('              Without it, root-counting by |xi_n|<1<=|xi_{n+1}| can falsely')
disp('              claim that the model is DET.')
disp('              This is what Chris Sims emphasizes in his gensys algorithm.')
disp('              Both this approach and gensys use the same QZ method')
disp('              based generalized Schur decomposition theorem. What was missing')
disp('              was to check the existence of the Omega_sn.')
disp('           2) Therefore, root-counting +checking the existence of the Omegs_sn')
disp('              becomes completely equivalent to gensys in classifying LRE models.')
disp('              "eu" checks existence of stable solution and its uniqueness.')
disp('              eu(1) =1 if a stable solution exists and eu(1)<0 otherwise.')
disp('              eu(2) =1 if a stable solution is unique and 0 if more than a stable solution exists.')
disp('              The conditions for DET/INDET/NSS precisely corresonds to eu=[1 1] / [1 0] / eu(1)<=0.')
disp('           3) Note that, However, Root-counting fails ONLY in the special case of DIA models')
disp('             in which Omega_ns does not exist. This can arise ONLY when ')
disp('             a model consists of completely decoupled equations.')
disp('  ')
disp('=========================================================')
disp('    QZ method now refers to as [2], which is idenitical to Root-Counting + Existence of Omega_ns.')
disp('    Now MOD method = QZ method = gensys')
disp('    Application of the MOD method therefore boils down to solving for the MOD solution.')
disp('    QZ method is a stand-alone procedure to fulfill the classification of LRE models because it automatically')
disp('    computes the MOD solution')
disp('    So the MOD method can be completed by QZ technique.')
disp('    The Groebner Basis (GB) technique proposed by FRWZ for MSRE models') 
disp('    plays exactly the same role as QZ method for LRE models, Theoretically')
disp('    ** The Difference is feasibility. GB identifies the MOD solution by')
disp('    computing ALL MSV solutions and find the most stable one.')
disp('    Unfortunately, the number of solutions increases exponentially ')
disp('    with the power commesurate with the dimension of the model.')
disp('    The computational time exponentially with the number of solutions.')
disp('    For this reason, it works only for a very small set of models.')
disp('    Therefore, there must be an efficient way of computing the MOD solution.')
disp('    MOD classification provides an extremely simple identification condition for the MOD solution in the case of')
disp('    which is the same as the condition for DA family.')
disp('To this end, we propose a sequential approach of applying the FM followed by GB if necessary.')
disp('   The original FM is known as picking up the MOD solution, except for a') 
disp('   very special models having a block-recursive structure with a special')
disp('   assumption. This is modified to pick up the MOD solution in Cho(2019).')
disp('   While there is no known proof for the equivalence between FS and MOD,')
disp('   there has been no exception so far for the DA models.')
disp('   The only cases in which the FS does not coincide with the MOD solution,')
disp('   are when the MOD solution is complex-valued, for which the FS cannot exist.')
disp('   This might sound like a bad news, but in fact it is a good news in the ')
disp('   sense that non-existence of the FS implies almost surely that the model')
disp('   has a complex-valued MOD solution which is by definition a DIA model, thus')
disp('   it cannot be determinate.')
disp('   Computation is  as comparable to standard solution method, once it exists.') 
disp('   The proposed methodology is nevertheless complete because one can always')
disp('   apply the GB approach.')

disp('  Of course, one does not need to use this sequential technique to find')
disp('  the MOD solution for LRE models because QZ method can always do it.')
disp('  However, it is not the case for MSRE models.')

disp('  The purpose of this code is to numerically show the equivalence of') 
disp('  This sequential approach with the QZ method and gensys algorithm.')
disp('  So that one can get a sense of feasibility of the sequential solution')
disp('  Technique.')


disp('========================================================')

disp('Three Methods for identifying the MOD solution and Classification')
disp('   [1] MOD Method')
disp('   [1-1] QZ method :  This is a stand-alone procedure using QZ decomposition.')
disp('                      Matlab Code: qzmlre.m ') 
disp('                      This corresponds to GB(Groebner Basis) Approach in MSRE model.')
disp('   [1-2] FM+QZ :      Using FM, followed by QZ if necessary.')
disp('                      This corresponds to modmsre_FM_GB Approach in MSRE model.')
disp('   [2] gensys Approach: This is a stand-alone procedure')
disp('                      Matlab Code: qzmlre2gensys.m : transform the model into gensys form and apply gensys.')

disp('  ')

disp('=============================================================================================')
disp(' This code tests the three different approachs using 8 Atheoretical Models') 
disp('   --Model--  -----Criteria---------  --------- FM First -----------   -------QZ or gensys----- ') 
disp('              DA?  Class?  real MOD?  FS Exist?  FS=MOD?   Class       QZ or gensys   Class  ')            
disp('                                                         identified?   required?    identified? ')         
disp('   Model=1111: DA   DET     Real    	Yes      Yes        Yes          No           Yes             ')            
disp('   Model=1112: DA   DET     Real    	Yes      Yes        Yes          No           Yes              ')          
disp('   Model=121 : DA   INDET   Real    	Yes      Yes        Yes          No           Yes              ')   
disp('   Model=131 : DA   NSS     Real    	Yes      Yes        Yes          No           Yes              ')   
disp('   Model=221 : DIA  INDET   Real     	Yes      unknown    Yes          No           Yes              ')     
disp('   Model=222 : DIA  INDET   Complex  	No       unknown    No           Yes          Yes              ')   
disp('   Model=231 : DIA  NSS     Real        Yes      unknown    No           Yes          Yes             ')    
disp('   Model=232 : DIA  NSS     Complex     No       unknown    No           Yes          Yes             ')    
disp('=============================================================================================')


disp('  ')
for CASE=1:8
    
if      CASE==1, Model=1111; disp('=== Model 1111 ( DA / DET / R ) ======================')
                 disp('  : Normal case: Not Block-Recursive, No decoupled equations')
                 disp('  ')
elseif  CASE==2, Model=1112; disp('=== Model 1112 ( DA / DET / R ) ======================') 
                 disp('  : Example 2 of Sims (2007) : Block recursive Model in which the explosive equation depends on the autonomous block.')
                 disp('  : Modified FM yields OmegaK = Omega_MOD. This was not in the Original FM.')
                 disp('  ')
elseif  CASE==3, Model=121; disp('=== Model 121 ( DA / INDET / R ) ======================')
                 disp('  : Model 1111 transformed to be indeterminate')
                 disp('  ')
elseif  CASE==4, Model=131; disp('=== Model 131 ( DA / NSS / R ) ======================')
                 disp('  : Model 1111 transformed to have NSS')
                 disp('  ')
elseif  CASE==5, Model=221; disp('=== Model 221 ( DIA / INDET / R ) ======================')
                 disp('  : Model 231 transformed to be INDET')
                 disp('  ')
elseif  CASE==6, Model=222; disp('=== Model 222 ( DIA / INDET / C ) ======================')
                 disp('  : Model 232 transformed to be INDET')
                 disp('  ')
                 
elseif  CASE==7, Model=231; disp('=== Model 231 ( DIA /NSS / R ) ======================') 
                 disp('  : Example 1 of Sims (2007) : Decoupled System with one explosive equation without forward-looking term')
                 disp('   and the other one is an autonomous block, which is an indeterminate model in isolation')
                 disp('  ')
elseif  CASE==8, Model=232; disp('=== Model 232 ( DIA /NSS / C ) ======================') 
                 disp('  : Example 3 of Sims (2007) : Omega_MOD is complex-valued, thus not unique.')    
                 disp('  ')
end

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
    n=size(A,1);
        
%% [MOD-A1]  
        [DETCqzm,Omega,Gamma,F,Geig,gvAll,OmegaAll,GammaAll,FAll]=qzmlre(A,B);
        

%% [MOD-A2]      
        % Step 1: Apply Modified Forward Method
            [DETC,FCC,OmegaK,GammaK,FK]=fmlre(A,B); 
        % Step 2: Apply QZ method if FM fails to detect DET/INDET/NSS    

%% [2] gensys Method
        % The following code transforms the model into gensys form and apply gensys.m
        [YG,YC,impact,fmat,fwt,ywt,gev,eu,loose]=qzmlre2gensys(A,B);   
        
%% [3] Root-Counting + existence of Omega(xi_1,...,xi_n): 
         
        

%% Results
disp(strjoin(['   [1-1] [r(OmegaMOD) r(FMOD)  r(OmegaMOD)*r(FMOD) ] by QZ =',string(DETCqzm)]))
             
                if      round(DETCqzm(3),7)<1 && DETCqzm(1)<1 && DETCqzm(2)<=1,      disp('        The model is DA / DET / R')
                elseif  round(DETCqzm(3),7)<1 && DETCqzm(1)<1 && DETCqzm(2)>1,       disp('        The model is DA / INDET / R')
                elseif  round(DETCqzm(3),7)<1 && DETCqzm(1)>1,                       disp('        The model is DA / INDET / R')
                elseif  round(DETCqzm(3),7)>=1 && DETCqzm(2)>1 && isreal(Omega),     disp('        The model is DIA / INDET / R')
                elseif  round(DETCqzm(3),7)>=1 && DETCqzm(2)>1 && ~isreal(Omega),    disp('        The model is DIA / INDET / C')
                elseif  round(DETCqzm(3),7)>=1 && DETCqzm(1)>=1 && isreal(Omega),    disp('        The model is DIA / NSS / R')
                elseif  round(DETCqzm(3),7)>=1 && DETCqzm(1)>=1 && ~isreal(Omega),   disp('        The model is DIA / NSS / C')
                end
if isnan(DETC(1))
    disp('   [1-2] MOD_FM_QZ: [r(OmegaK) r(FK)  r(OmegaK)*r(FK) ] by FM = NaN NaN NaN')
else
    disp(strjoin(['   [1-2] MOD_FM_QZ: [r(OmegaK) r(FK)  r(OmegaK)*r(FK) ] by FM =',string(DETC)]))
end

                if      round(DETC(3),7)<1 && DETC(1)<1 && DETC(2)<=1,          disp('        r(OmegaK)*r(FK)<1, thus OmegaK=Omega_MOD. The model is DA / DET / R')
                elseif  round(DETC(3),7)<1 && DETC(1)<1 && DETC(2)>1,           disp('        r(OmegaK)*r(FK)<1, thus OmegaK=Omega_MOD. The model is DA / INDET / R')
                elseif  round(DETC(3),7)<1 && DETC(1)>1,                        disp('        r(OmegaK)*r(FK)<1, thus OmegaK=Omega_MOD. The model is DA / NSS / R')
                elseif  round(DETC(3),7)>=1 && DETC(2)>1 && isreal(OmegaK),      disp('        OmegaK=OmegaMOD is not verifiable. But The model is ?? / INDET /??') 
                                                                                disp('        Not needed to apply QZ. But applying it confirms that OmegaK=OmegaMOD as')
                                                                                disp(strjoin(['        [r(OmegaMOD) r(FMOD)  r(OmegaMOD)*r(FMOD) ]=',string(DETCqzm)]))
                                                                                disp('        Therefore, the model is DIA / INDET / R')
                                                                                r1MOD=max(abs(eig(Omega))); 
                                                                                disp(strjoin(['        Omega_ns does NOT exist because r(OmegaMOD)=',string(r1MOD),'>|xi_n|=',string(Geig(n))])) 
                                                                                disp('        This shows the importance of checking the existence of Omega_ns, which is done by QZ method.')
               elseif  round(DETC(3),7)>=1 && DETC(1)>=1 && isreal(OmegaK),     disp('        OmegaK=OmegaMOD is not verifiable. Since r(OmegaK)>1, apply QZ method.') 
                                                                                disp(strjoin(['        [r(OmegaMOD) r(FMOD)  r(OmegaMOD)*r(FMOD) ]=',string(DETCqzm)]))
                                                                                disp('        OmegaK=OmegaMOD is confirmed and thus, the model is DIA / INDET / R')
                                                                                disp('        Therefore, the model is DIA / NSS / R')
                                                                                disp('        This shows the importance of checking the existence of Omega_ns, which is done by QZ method.')  
                elseif  isnan(OmegaK)
                    if round(DETCqzm(3),7)>=1 && DETCqzm(2)>1 && ~isreal(Omega)
                                                                                disp('        OmegaK does not exist. Apply QZ method.')
                                                                                disp(strjoin(['        [r(OmegaMOD) r(FMOD)  r(OmegaMOD)*r(FMOD) ]=',string(DETCqzm)]))
                                                                                disp('        The MOD solution OmegaMOD='), disp(Omega)
                                                                                disp('        The MOD solution is complex-valued and r2(MOD)>1, thus the model is DIA / INDET / C.')
                    
                    elseif  round(DETCqzm(3),7)>=1 && DETCqzm(1)>=1 && ~isreal(Omega)   
                                                                                disp('        OmegaK does not exist. Apply QZ method.')
                                                                                disp(strjoin(['        [r(OmegaMOD) r(FMOD)  r(OmegaMOD)*r(FMOD) ]=',string(DETCqzm)]))
                                                                                disp('        The MOD solution OmegaMOD='), disp(Omega)
                                                                                disp('        The MOD solution is complex-valued and r1(MOD)>=1. The model is therefore, DIA / INDET / C.')
                    end
                end
disp('           ')

%% gensys 
                if eu==[1;1],                       disp(strjoin(['   [2] gensys : The model is DET because eu=',string(eu')]))
                elseif eu==[1;0],                   disp(strjoin(['   [2] gensys : The model is INDET because eu=',string(eu')]))
                elseif eu(1)<=0,                    disp(strjoin(['   [2] gensys : The model has NSS because because eu=',string(eu')]))
                end

end