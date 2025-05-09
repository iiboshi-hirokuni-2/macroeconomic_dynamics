%
%  Chapter 3
%  Sec 3.4.1
%

clc
clear all

% set parameter
P.beta = 0.98;
P.theta = 0.36 ;
P.delta = 0.1;
P.A = 0.5;
P.tau =0;
tau_table = [ 0, 0.1, 0.2];

% initial value
K = 3.5;
C = 0.8;
H = 0.6;
x0= [K,C,H];

% solve steady state
R0 = SS(x0,P);

disp( ' '); disp( ' '); disp( ' ');
disp('Table 3.1')

for i = 1:3
     P.tau=tau_table(i);     
     
     % options = optimoptions('fsolve','Display','iter');
     options = optimoptions('fsolve','Display','off');
     fun = @(x)SS(x,P);
     [x, f ]= fsolve(fun,x0,options);
     c = x(2);
     h = x(3);

     x= real(x);
     u = log(c)+P.A*log(1-h);

     disp( ' ');
     disp([ 'tau = ', num2str(P.tau) ]);
     disp( [ 'C_bar=', num2str(x(2))  ]);
     disp( [ 'K_bar=', num2str(x(1))  ]);     
     disp( [ 'H_bar=', num2str(x(3))  ]);
     disp( [ 'u_bar = ', num2str(real(u))]);

end     

%% function

function [ R ] = SS(x, P)
   R= zeros(3,1); % residual 

   K=x(1); C=x(2); H=x(3); 
   
   % page 48
   Y = K^P.theta*H^(1-P.theta); % Production Function

   % FOC 1
   R(1) = P.delta*K - ( Y - C );
   
   % FOC 2
   R(2) = P.beta*(P.theta*K^(P.theta-1)*H^(1-P.theta)+(1-P.delta)) -1;

   % FOC 3
   R(3)= -(-P.A/(1-H))/(1/C) - (1-P.tau) * K^P.theta*(1-P.theta)*H^(-P.theta);

   % disp( num2str( (-P.A/(1-H))/(1/C) ) );
   % disp(num2str(-(1-P.tau) * K^P.theta*(1-P.theta)*H^(-P.theta)));
end
