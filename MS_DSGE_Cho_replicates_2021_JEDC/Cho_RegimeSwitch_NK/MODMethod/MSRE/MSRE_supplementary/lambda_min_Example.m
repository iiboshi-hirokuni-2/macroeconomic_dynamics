%% Testing lambda_min
%% Analytical Solution for Lambda
disp(' =================== Problem ===========================')
disp(' For given P and F{s(t),s(t+1)} such that xi2=r(Phi_kron(F,F)),')
disp(' Proposition 1 of Cho(2018) shows that there exists a Lambda{s(t),s(t+1)}') 
disp(' such that tau2=r(BarPhi_kron(Lambda,Lambda))=1/xi2 ')
disp(' That is, tau2*xi2=1')
disp(' This code construct an arbitrary F, and find such a Lambda.')
disp('   ')

clear

% (1) Size of the Model and Number of Regime

    n=4; % Dimension of F(s(t),s(t+1))
    S=2; % # of regimes
    
    % Construct a random transition probability matrix 
    P0=rand(S,S);
        %Arbitrary zeros
            %P0(2,1)=0; 
        %Arbitrary Absorbing state 
            %P0(1,2:end)=0;   
        for i=1:S, P(i,:)=P0(i,:)/sum(P0(i,:)); end 
    
    % Construct a randome F with rank<n in some states
        for i=1:S, k{i,1}=randi(n); end
        %for i=1:S, k{i,1}=n; end
        
        for i=1:S,  SeedV{i}=orth(randn(n,k{i}));  end
        for i=1:S, for j=1:S
            F{i,j}=SeedV{i}*SeedV{i}'*rand(n,n); 
        end, end 
        % Assign zeros to states (s(t),s(t+1))=(i,j) if P(i,j)=0; 
        for i=1:S, for j=1:S
              if P(i,j)==0, F{i,j}=zeros(n,n); end
        end, end 
  
    % Given P and F, find Lambda_min
    [Lambda_min,tau2min,xi2min,V,Phi,D,u,Q]=lambda_min(P,F);

        disp('[xi2 tau2, xi2*tau2] is')
        disp([xi2min tau2min xi2min*tau2min])
        disp('Type F and Lambda_min to see the results.')





    
    