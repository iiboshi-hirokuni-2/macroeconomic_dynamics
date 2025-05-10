function [DETC,FCC,OmegaK,FK,DETCMOD,OmegaMOD,FMOD]=modlre_FM_QZ(A,B)
%% This is a general code implementing the MOD method in the most
%% general way using the forward method(FM), and the QZ method for LRE model if necessary.
%% Acronyms
% DA / DIA : Determinacy-Adimissible / Determinacy-Inadmissible
% DET/INDET/NSS : Determinate / Indeterminate / No Stable Solution
% FM / QZ : Forward Method / QZ method
% FS / MOD : Forward Solution / MOD solution
%% This is the LRE analogue to modmsre_FM_GB for MSRE models.
%   This code completes the MOD method. 
%   Unlike the MSRE models, QZ method (qzmlre.m) is extremely efficient.
%   Therefore, this code and using qzmlre.m alone are equivalent. 
%   These are also equivalent to gensys algorithm, but reformulated in our 
%   model and solution format.
%   Note however, it is not the case for MSRE models. GB approach plays the 
%   same role as QZ method, but extremely inefficient. Therefore, the
%   hybrid approach of modmsre_FM_GB.m is a must. 
%   This code can be used to show the intuition behind the hybrid approach.
      
%% Algorithm: FM, followed by QZ method if necessary  
    Opt.maxK=1000;
    [DETC,FCC,OmegaK,~,FK]=fmlre(A,B,[],[],Opt);
        disp('========================================================')
        disp('The following acronyms implies that ')
        disp('DETC contains the following information for a given solution Omega^{i}(s(t))')
        disp('    DETC(1)=r(Omega^{i},Omega^{i})')
        disp('    DETC(2)=r(F^{i},F^{i})')
        disp('    DETC(3)=DETC(1)*DETC(2)')
        disp('DA=Determinacy-Admissible, DIA=Determinacy-Inadmissible')
        disp('Model is DA if DETC(3)<1 for the MOD solution.')
        disp('    ')
        disp('========================================================')
        disp('Results from Forward Method (FM) =======================')
        disp('========================================================')  
   	if ~isnan(DETC(3))
        disp('    DETC for the forward solution (FS) OmegaK'), disp(DETC)
        disp('    The forward solution =')
        disp(OmegaK)
    end

    if DETC(3)<1
        disp('    DETC(3)<1, thus FS = MOD solution and the model is DA.')
        disp('    No need for QZ method')     
        DETCMOD=DETC;OmegaMOD=OmegaK; FMOD=FK;
        return
    elseif DETC(3)>=1 && DETC(1)<1
        disp('    DETC(3)>=1, thus FS may or may not be the MOD solution.')
        disp('    Nevertheless, DETC(1)<1 and DETC(2)>1 implies that')
        disp('    the model is INDET regardless of FS=MOD solution.')
        disp('    To see QZ results anyway, use qzmlre(A,B) instead of this code.')
        DETCMOD ="Not confirmed"; OmegaMOD="Not confirmed"; FMOD="Not confirmed";
        return
    elseif DETC(3)>=1 && DETC(1)>=1 
        disp('    DETC(3)>=1, thus FS may or may not be the MOD solution.')  
        disp('    Applying the QZ method to identify MOD solution.')
        [DETCMOD,OmegaMOD,~,FMOD]=qzmlre(A,B); 
        disp('========================================================')
        disp('Results from QZ Method =======================')
        disp('========================================================')  
        disp('    DETMOD from QZ method'), disp(DETCMOD)
        disp('    The MOD solution =')
        disp(OmegaMOD)
    elseif isnan(DETC(3))
        disp('    FS does not exist.')
        disp('    You may set a higher Opt.maxK at the beginning of modlre_FM_QZ.m.')
        disp('    Applying the QZ method to identify MOD solution.')
        [DETCMOD,OmegaMOD,~,FMOD]=qzmlre(A,B); 
        disp('========================================================')
        disp('Results from QZ Method =================================')
        disp('========================================================')  
        disp('    DETMOD from QZ method'), disp(DETCMOD)
        disp('    The MOD solution =')
        disp(OmegaMOD)
        
    end
    
