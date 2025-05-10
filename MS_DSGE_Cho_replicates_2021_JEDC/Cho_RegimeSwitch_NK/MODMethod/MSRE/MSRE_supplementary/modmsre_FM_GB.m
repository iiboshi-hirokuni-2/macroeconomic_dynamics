function [DETC,FCC,OmegaK,FK,DETCMOD,OmegaMOD,FMOD,DETC_All,AllOmegas]=modmsre_FM_GB(P,A,B) 
%% This is a general code implementing the MOD method in the most
%% general way using the forward method(FM), and the Groebner basis(GB) approach if necessary
%% Important: To use this code, users must install Singular program.
%% Refer to gbmsre.m for GB for more detail.
%% Tradeoff
%   (1) gbmsre.m 
%       Pros: This is a complete stand-alone code for implementing the MOD method.
%       Cons: It works only for very small scale models only due to the
%               computational burden.
%   (2) fmmsre.m 
%       Pros: - This works almost surely for the class of DA (Determinacy-Admissible) models and some
%               indeterminate DIA (Determinacy-Inadmissible) models.
%             - Computationally fast engough, comarable to standard methods for LRE models.
%       Cons: - For some DIA models, it does not work.
%             - It is not proved theoretically that the resulting forward solution = MOD
%               solution for DA models.
%% Evaluation
%   - Efficient implementation of the MOD method is to first try FM
%     followed by GB if FM fails to identify MOD solution.
%   - Virtually all economic models belong to DA models.
%   - Forward solution = MOD solution for Virtually all models: Not a
%     single counterexample has been found for DA models.
%% Usage 
%   (1) Install Singular
%   (2) Write your n-dimensional MSRE model with S regimes in terms of P, A and B and run it.
%   (3) Follow the instruction displayed in matlab command.

%% Algorithm: FM, followed by GB if necessary  
    Opt.maxK=1000;
    [DETC,FCC,OmegaK,~,FK]=fmmsre(P,A,B,[],[],Opt);
        disp('========================================================')
        disp('The following acronyms implies that=====================')
        disp('DETC contains the following information for a given solution Omega^{i}(s(t))')
        disp('    DETC(1)=r(BarPsi(kron(Omega^{i},Omega^{i})))')
        disp('    DETC(2)=r(Psi(kron(F^{i},F^{i})))')
        disp('    DETC(3)=DETC(1)*DETC(2)')
        disp('    DETC(4)=r(Psi(kron((Omega^{i})^T,F^{i})))')
        disp('DA=Determinacy-Admissible, DIA=Determinacy-Inadmissible')
        disp('If DETC(3)<1 and real-valued, the solution is the MOD solution and the Model is DA.')
        disp('    ')
        disp('========================================================')
        disp('Results from Forward Method (FM) =======================')
        disp('========================================================')  
   	if ~isnan(DETC(3))
        disp('    DETC for the forward solution, OmegaK(s(t))'), disp(DETC)
        disp('    The forward solution at each regime s(t)=')
        for s=1:size(P,1), disp(OmegaK{s,1}), end
    end

    if DETC(3)<1
        disp('    DETC(3)<1, thus FS = MOD solution and the model is DA.')
        disp('    No need for GB Approach')     
        DETCMOD=DETC;OmegaMOD=OmegaK; FMOD=FK;
        DETC_All="Not Available"; AllOmegas="Not Available"; 
        return
    elseif DETC(3)>=1 && DETC(1)<1
        disp('    DETC(3)>=1, thus FS may or may not be the MOD solution.')
        disp('    Nevertheless, DETC(1)<1 and DETC(2)>1 implies that')
        disp('    the model is INDET regardless of FS=MOD solution.')
        disp('    To see GB results anyway, use gbmsre(P,A,B) instead of this code.')
        DETCMOD ="Not confirmed"; OmegaMOD="Not confirmed"; FMOD="Not confirmed";
        DETC_All="Not Available"; AllOmegas="Not Available"; 
        return
    elseif DETC(3)>=1 && DETC(1)>=1 
        disp('    DETC(3)>=1, thus FS may or may not be the MOD solution.')  
        disp('    Applying the GB approach to identify MOD solution.')
        disp('    Follow the instruction below')
        [DETCMOD,OmegaMOD,FMOD,DETC_All,AllOmegas]=gbmsre(P,A,B); 
        disp('    Since DETCMOD(3)>=1, the model is DIA.')
        disp('    INDET/NSS can be identified by DETCMOD(1) and DETCMOD(2).')
    elseif isnan(DETC(3))
        disp('    FS does not exist.')
        disp('    You may set a higher Opt.maxK at the beginning of modmsre_FM_GB.m.')
        disp('    Applying the GB approach to identify MOD solution.')
        disp('    Follow the instruction below')
        [DETCMOD,OmegaMOD,FMOD,DETC_All,AllOmegas]=gbmsre(P,A,B); 
        disp('    Since the MOD solution is complex-valued, the model is DIA.')
        disp('    INDET/NSS can be identified by DETCMOD(1).')
    end

    