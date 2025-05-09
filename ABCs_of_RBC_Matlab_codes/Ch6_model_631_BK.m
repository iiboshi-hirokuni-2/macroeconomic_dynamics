
clc
clear all

tic

%Define parameters
syms THETA BETTA DELTA AA GAMMA
syms  r_ss h_ss y_ss c_ss k_ss

%Define variables 
syms ct  rt  yyt  ht  kp lamt;  % control variables
syms cf rf yyf hf kpf lamf;     % future period 
syms cb rb yyb hb kt  lamb;     % state and lagged variables   
syms eps;

%% Set FOC of Model: Page 100 
eq1 = ct - cf + BETTA*r_ss*rf ;

eq2 = yyt - ht/(1-h_ss)-ct;

eq3 = y_ss*yyt - c_ss*ct + k_ss*((1-DELTA)*kt - kp);

eq4 = lamt + THETA*kt + (1-THETA)*ht - yyt ;

eq5 = yyt - kt - rt;

eq6 = -lamt + GAMMA*lamb + eps; 

%% Create function f
f = [eq1;eq2; eq3; eq4; eq5; eq6];

xf = [kp, lamt, yyt, ht, cf, rf]; 
x =  [kt, lamb, yyb hb, ct, rt]; 
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
AA = 1.72;
para(1)=BETTA; para(2)=THETA; para(3)=DELTA; para(4)= AA;

%% solve steatdystate
%x=  [kp, y, c h, r];
x0= [0, 0, 0, 0, 0,0];
fun = @(x)steatdystate(x,para);
x1= fsolve(fun,x0);

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
eigen4 = diag(S(1:5,1:5)./T(1:5,1:5))

% unstable  page 135
Z= Z';
N=inv(Z(5:6,5:6))*Z(5:6,1:4)

% stable 
B11= B(1:4,1:4);
B12=B(1:4,5:6);
A11 = A(1:4,1:4);
A12 = A(1:4,5:6);

(B11-B12*N)
(A11-A12*N)  

 RR  =  inv(B11-B12*N)*(A11-A12*N)

%% Impulse response
 % x =  [kt, lamb, yyb hb, ct, rt]; 

 eps = 0.01;
 time = 100;
 
 x=zeros(6,time);
 x(:,1)=[0; eps; 0;0;0;0];
   
 for  t =2:time 
       x(1:4,t) =  RR * x(1:4,t-1);
       x(5:6,t) =  -1*N * x(1:4,t-1);

 end

%%
% y =[ 'K','lamdda','Y','H','C','R']

figure('Name','IRF')
l1=plot(x(2,:)','b-','Linewidth',2);
hold on  
   l2=plot(x(1,:)','r-','Linewidth',2);
   l3=plot(x(3,:)','g-.','Linewidth',2);
   l4=plot(x(5,:)','m-.','Linewidth',2);
   l5=plot(x(4,:)','c-','Linewidth',2);
   l6=plot(x(6,:)','k:','Linewidth',2);
hold off
set(gca, 'Fontsize',12)
legend([l1, l2, l3, l4, l5,l6],...
      {'lamdda','k','y', 'c' ,'h','r'},'Fontsize',14)
title('Impulse response (BK)','Fontsize',16)

disp([ 'cal time =' num2str(toc) 'sec' ])

%%
function R = steatdystate(x,para)

BETTA=para(1); THETA=para(2); DELTA=para(3); A=para(4);

% [kp, y, c h, r lam];
k=x(1); y=x(2); c=x(3); h=x(4); r =x(5); lam=x(6);
cp=x(3); rp=x(5); kp=x(1);


R(1) = BETTA*(exp(c)/exp(cp))*(exp(rp)+(1-DELTA)) -1 ;

R(2) = (1- THETA)*(1-exp(h))*exp(y)/exp(h)- A*exp(c);

R(3) = exp(y) + (1-DELTA)*exp(k) - exp(kp)- exp(c);

R(4) = exp(lam)*exp(k)^THETA*exp(h)^(1-THETA) -exp(y) ;

R(5) = THETA*exp(y)/exp(k) - exp(r);

end

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
