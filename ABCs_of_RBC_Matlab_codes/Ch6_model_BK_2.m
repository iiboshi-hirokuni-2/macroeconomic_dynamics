
 clc

tic

%Define parameters
syms THETA BETTA DELTA A GAMMA;
syms  r_ss h_ss y_ss c_ss k_ss ;

%Define variables 
syms c  r  y  h  kp lam;  % control variables
syms cf rf yf hf kpf lamf; % future period 
syms cb rb yb hb k  lamb; % state and lagged variables   
syms eps

%% Set FOC of Model: Page 100 
eq1= y_ss*y - c_ss*c + k_ss*((1-DELTA)*k - kp);

 % eq2 = y - h/(1-h_ss)-c;
 eq2 = lam - GAMMA *lamb - eps;

 % eq4 = lam + THETA*k + (1-THETA)*h - y ;
 eq3 = lam - THETA*y + THETA*k - (1-THETA)*c  ;

eq4 = -y + k + r;

eq5 = -c + cf - BETTA*r_ss*rf ;

%% Create function f
f = [eq1;eq2; eq3; eq4; eq5];

xf = [kp, lam, y, cf, rf]; 
x =  [k, lamb, yb c, r]; 
z=    eps; % exogenous variables

n = size(f,1); % number of equations
nx = size(x,1); % number of variables

%% Compute the first  derivatives of f
fxf=jacobian(f,xf);
fx=jacobian(f,x);
fz=jacobian(f,z);


% set parameters
BETTA = 0.99;
THETA = 0.36;
DELTA = 0.025;
GAMMA = 0.95;
A = 1.72;
para(1)=BETTA; para(2)=THETA; para(3)=DELTA; para(4)= A;

%% solve steatdystate
%x=  [kp, y, c h, r];
x0= [0, 0, 0, 0, 0,0];
fun = @(x)steatdystate(x,para);
x1= fsolve(fun,x0);
x1= real(x1);

k_ss=exp(x1(1)); y_ss=exp(x1(2)); c_ss=exp(x1(3)); h_ss=exp(x1(4)); r_ss =exp(x1(5));

disp('Steady State :')
disp([ 'k_ss= ', num2str(exp(x1(1)))]);
disp([ 'y_ss= ', num2str(exp(x1(2)))]);
disp([ 'c_ss= ', num2str(exp(x1(3)))]);
disp([ 'h_ss= ', num2str(exp(x1(4)))]);
disp([ 'r_ss= ', num2str(exp(x1(5)))]);
disp([ 'lam_ss= ', num2str(exp(x1(6)))]);

disp(' '); 

%% 

disp('df/dxf (B) =')
disp(fxf)
num_fxf = eval(fxf);  B=-1*num_fxf;
disp('df/dxf (B) =');
disp(num_fxf);

disp('df/dx (A) =')
disp(fx)
num_fx = eval(fx);  A=num_fx;
disp('df/dx (A) =');
disp(num_fx);


% inverse matrix of B 
inv(B)


% [N,L,C,D,alphabeta]=modelschur(A,B)


%%  P139:  QZ  Decomposition  Sec 6.8.4

[T,S,q,z] = qz(B,A);   % upper triangular factorization of the matrix pencil b-za
     stake=1;
[T,S,q,Z] = qzdiv(stake,T,S,q,z);   % reordering of generalized eigenvalues in ascending order

eigen5 = diag(S./T)
eigen4 = diag(S(1:4,1:4)./T(1:4,1:4))

% unstable  page 135
Z= Z';
N=inv(Z(4:5,4:5))*Z(4:5,1:3)

% stable 
B11= B(1:3,1:3);
B12=B(1:3,4:5);
A11 = A(1:3,1:3);
A12 = A(1:3,4:5);

(B11-B12*N)
(A11-A12*N)  

 RR  =  inv(B11-B12*N)*(A11-A12*N)

 %% Impulse response
 
 eps = 0.01;
 time = 100;
 
 x=zeros(5,time);
 x(:,1)=[0; eps; 0;0;0];
   
 for  t =2:time 
       x(1:3,t) =  RR * x(1:3,t-1);
       x(4:5,t) =  -1*N * x(1:3,t-1);

 end

figure(200)
l1=plot(x(1,:)','g--','Linewidth',2);
hold on  
   l2=plot(x(2,:)','m-','Linewidth',2);
   l3=plot(x(3,:)','r-.','Linewidth',2);
   l4=plot(x(4,:)','b-.','Linewidth',2);
   l5=plot(x(5,:)','c-','Linewidth',2);
hold off
legend([l1, l2, l3, l4, l5],...
       {'K','lamdda','Y','C','R' })
title('Impulse response (BK)','Fontsize',16)

xlabel('Time')
set(gca, 'Fontsize',12)

disp([ 'cal time =' num2str(toc) 'sec' ])


%%
%% function
function [A,B,Q,Z] = qzdiv(stake,A,B,Q,Z)
%
% Written by Chris Sims
%
% Takes U.T. matrices A, B, orthonormal matrices Q,Z, rearranges them
% so that all cases of abs(B(i,i)/A(i,i))>stake are in lower right 
% corner, while preserving U.T. and orthonormal properties and Q'AZ' and
% Q'BZ'.
%
[n jnk] = size(A);
root = abs([diag(A) diag(B)]);
root(:,1) = root(:,1)-(root(:,1)<1.e-13).*(root(:,1)+root(:,2));
root(:,2) = root(:,2)./root(:,1);
for i = n:-1:1
   m=0;
   for j=i:-1:1
      if (root(j,2) > stake | root(j,2) < -.1) 
         m=j;
         break
      end
   end
   if (m==0) 
      return 
   end
   for k=m:1:i-1
      [A B Q Z] = qzswitch(k,A,B,Q,Z);
      tmp = root(k,2);
      root(k,2) = root(k+1,2);
      root(k+1,2) = tmp;
   end
end   

end % end of function
%==========================================%
function [A,B,Q,Z] = qzswitch(i,A,B,Q,Z)
%function [A,B,Q,Z] = qzswitch(i,A,B,Q,Z)
% Written by Chris Sims
% Takes U.T. matrices A, B, orthonormal matrices Q,Z, interchanges
% diagonal elements i and i+1 of both A and B, while maintaining 
% Q'AZ' and Q'BZ' unchanged.  Does nothing if ratios of diagonal elements
% in A and B at i and i+1 are the same.  Aborts if diagonal elements of
% both A and B are zero at either position.
%

a = A(i,i); d = B(i,i); b = A(i,i+1); e = B(i,i+1);
c = A(i+1,i+1); f = B(i+1,i+1); 
wz = [c*e-f*b, (c*d-f*a)'];
xy = [(b*d-e*a)', (c*d-f*a)'];
n = sqrt(wz*wz');
m = sqrt(xy*xy');

if n == 0
   return
else
   wz = n\wz;
   xy = m\xy;
   wz = [wz; -wz(2)', wz(1)'];
   xy = [xy;-xy(2)', xy(1)'];
   A(i:i+1,:) = xy*A(i:i+1,:);
   B(i:i+1,:) = xy*B(i:i+1,:);
   A(:,i:i+1) = A(:,i:i+1)*wz;
   B(:,i:i+1) = B(:,i:i+1)*wz;
   Z(:,i:i+1) = Z(:,i:i+1)*wz;
   Q(i:i+1,:) = xy*Q(i:i+1,:);
end

end % end of function





%%
function R = steatdystate(x,para)

BETTA=para(1); THETA=para(2); DELTA=para(3); A=para(4);

% [kp, y, c h, r lam];
k=x(1); y=x(2); c=x(3); h=x(4); r =x(5); lam=1;
cp=x(3); rp=x(5); kp=x(1);

R(1) = exp(h)-1/(1+A/(1-THETA)*(1-DELTA*THETA*BETTA/(1-BETTA*(1-DELTA)))) ;

R(2) = exp(k)-(THETA*BETTA/(1-BETTA*(1-DELTA)))^(1/(1-THETA))*exp(h);

R(3) = exp(y) - DELTA*exp(k) - exp(c);

R(4) = exp(lam)*exp(k)^THETA*exp(h)^(1-THETA) -exp(y) ;

R(5) = THETA*exp(lam)*exp(k)^(THETA-1)*h^(1-THETA) - exp(r);

end

