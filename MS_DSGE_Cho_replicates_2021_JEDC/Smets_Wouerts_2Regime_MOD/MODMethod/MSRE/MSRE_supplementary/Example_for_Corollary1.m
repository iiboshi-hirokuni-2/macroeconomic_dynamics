%% This is an numerical example for Proposition 4 of Cho(2019)
%% Uniqueness of Mean-square stable (MSS) solution does NOT imply Determinacy.
%% Whenever models contain lagged engogenous variables so there are
%% multiple MSV solutions.
%% Intuition 
%   Consider any Determinacy-Admissible(DA) model. 
%   Let Omega1, Omega2 be the MOD and the second most stable solutions in the MSS sense.
%   Proposition 4 comes from the fact that DETC1(2)*DETC2(1)>1 in general in MSRE models.
%   The problem is that this is true for DA models, not just DIA models.
%   In contrast, it is always the case that DETC1(2)*DETC2(1)=1 in LRE models as long as
%   no completely decoupled equations with special structures are included.
%   (Such a technical possibility is examined in Model 231, which is a DIA
%   model.) 
%   So consider any MSRE model with DETC2(1)*DETC1(2)>1. Then,
%   One can always make the model has a unique stable MSV solution,
%   but indeterminate by setting Atil( )=A(  )/alpha, Btil(  )=alpha*B(  )
%   where alpha=[Psi_kron(F^MOD,F^MOD)/Barpsi_kron(Omega^(2),Omega^(2))]^(1/4).
%   This transformation applies to both DA and DIA models.
%   For DIA model, Proposition is obvious, so this example shows that it
%   holds for DA model as well.

clear
%% Conisder any DA model. It can be DET/INDET/NSS.
%% The folowing trasformation makes it have a unique MSS solution but indeterminate.
    n=2; S=2;  
    P11=0.95; P22=0.9;
    P=[P11  1-P11;1-P22 P22];
    for k=1:1000
        for Si=1:S, B{Si,1}=randn(n,n); 
                    if n>1, B{Si,1}(:,1:end-1)=zeros(n,n-1); end
                    % To simplify the arbitrary model as GB cannot compute
                    % solutions if n>2 in general. 
            for Sj=1:S, A{Si,Sj}=randn(n,n); end
        end
        DETC=fmmsre(P,A,B);
        if DETC(3)<1 && DETC(2)>1, break, end % This is to find any DA model.
    end
    
    disp('==============================================================')
    disp('This code shows that uniqueness of MSV solution does not ')
    disp('imply determinacy in MSRE models.')
    disp('==============================================================')
    disp('First, an arbitrary Determinacy-Admissible model is generated. DETC is')
    disp(DETC);
    
        [DETCMOD,OmegaMOD,FMOD,DETC_All,AllOmegas,DETC_OmjtkF1]=gbmsre(P,A,B); 
        
    if DETC_All(2,1)>1 && DETC_All(1,2)>1, return, end    
 
        
    % Constructing an example such that uniqueness of stable MOD does not mean determinacy.    
     alpha=(DETC_All(1,2)/DETC_All(2,1))^(1/4);
     r12=DETC_All(2,1)*DETC_All(1,2);
     r12str="Barpsi_kron(Omega2,Omega2)*Psi_kron(F1,F1)=";
     disp(strjoin([r12str, string(r12)]))
     if round(r12,7)<=1, 
         disp(strjoin(['Since', r12str, '<=1', 'run this code again to try another model']))
         return
     else
        disp(strjoin(['Since', r12str, '>1', 'the following transformed model']))
        disp('has a unique stable MOD solution, but the model is indeterminate.')
        disp('The model is Atil(  )=A( )/alpha and Btil=B(  )*alpha where alpha=')
        disp(alpha)
     end
     
%% Result for the transformed model  
    DETCMOD_tanal=DETCMOD; 
    DETCMOD_tanal(1,1)=DETCMOD(1,1)*alpha^2; 
    DETCMOD_tanal(1,2)=DETCMOD(1,2)/alpha^2;
    DETC_All_tanal=DETC_All;
    DETC_All_tanal(:,1)=DETC_All(:,1)*alpha^2;
    DETC_All_tanal(:,2)=DETC_All(:,2)/alpha^2;
    ns=size(AllOmegas,1);
    for Si=1:S, OmegaMOD_tanal{Si,1}=OmegaMOD{Si,1}*alpha;
             for nsi=1:ns
                 AllOmegas_tanal{nsi,Si}=AllOmegas{nsi,Si}*alpha;   
             end
        for Sj=1:S, FMOD_tanal{Si,Sj}=FMOD{Si,Sj}/alpha; end
    end 
    %DETCMOD_tanal,OmegaMOD_tanal,FMOD_tanal,DETC_All_tanal,AllOmegas_tanal

        Sol_order=string((1:1:ns)');
        Sol_order(1)="MOD"; 
        if isequal(DETC_All_tanal(1,1),DETC_All_tanal(2,1)), Sol_order(2)="MOD"; end

        mdi=1; % Do not display if DETC(3)>>1000.
        for i=1:ns,
            if max(DETC_All_tanal(i,:))>=1000, mdi=i-1; break, 
            else, mdi=ns; 
            end 
        end

        D5str=string(DETC_All_tanal(:,5)); 
        for nsi=1:ns
            if DETC_All_tanal(nsi,5)==1, D5str(nsi,1)='real_valued'; 
            elseif DETC_All_tanal(nsi,5)==0, D5str(nsi,1)='complex_valued';
            end
        end
        disp('    Table below shows DETC(1) through DETC(4) information and real-valuedness')  
        disp('    for the transformed model.')
        disp('    ')
        varNames = {'Solution_Order','DETC','Solution_is'};
        T2=table(Sol_order(1:mdi,1),DETC_All_tanal(1:mdi,1:4),D5str(1:mdi,1),'VariableNames',varNames);
        disp(T2)
        disp('    The MOD solution  Omega(s) is')
        for s=1:S, disp(AllOmegas_tanal{1,s}), end
    
%% The result can also be checked by applying GB using the transformed model.
%    for Si=1:S, Btil{Si,1}=alpha*B{Si,1};
%        for Sj=1:S, Atil{Si,Sj}=A{Si,Sj}/alpha; end
%    end
%         [DETCMODtil,OmegaMODtil,FMODtil,DETC_Alltil,AllOmegastil]=gbmsre(P,Atil,Btil);   
%disp('The original model is determinacy-admissible.')
%disp('The transformed model has a unique stable MSV (MOD) solution')
%disp('But the transformed model is indeteminate as D1(MOD)<1 and D2(MOD)>1')

 