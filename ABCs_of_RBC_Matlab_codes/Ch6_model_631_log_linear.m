    
clc
clear all

tic

%Define parameters
syms THETA BETTA DELTA AA
syms  r_ss h_ss y_ss c_ss k_ss

%Define variables 
syms ct  rt  yyt  ht  kp lamt;  % control variables
syms cf rf yyf hf kpf lamf;     % future period 
syms cb rb yyb hb kt  lamb;     % state and lagged variables   


%% Set FOC of Model: Page 100 
eq1 = ct - cf + BETTA*r_ss*rf ;

eq2 = yyt - ht/(1-h_ss)-ct;

eq3 = y_ss*yyt - c_ss*ct + k_ss*((1-DELTA)*kt - kp);

eq4 = lamt + THETA*kt + (1-THETA)*ht - yyt ;

eq5 = yyt - kt - rt;

%% Create function f
f = [eq2; eq3; eq4; eq5];

% State Variables
xf = [kpf];
x = [kp]; 
xb = [kt];

% control variables
yf = [yyf, cf, hf, rf];
yt = [ yyt, ct ht, rt];
yb= [ yyb, cb hb, rb]; 

% exogenous variables
z= lamt; 
zf = lamf;

n = size(f,1); % number of equations
nx = size(x,1); % number of variables

%% Compute the first  derivatives of f
fx=jacobian(f,x);
fxb=jacobian(f,xb);
fyt=jacobian(f,yt);
fz=jacobian(f,z);

eq1xf=jacobian(eq1,xf);
eq1x =jacobian(eq1,x);
eq1xb =jacobian(eq1,xb);
eq1yf=jacobian(eq1,yf);
eq1yt =jacobian(eq1,yt);
eq1zf =jacobian(eq1,zf);

% set parameters
BETTA = 0.99;
THETA = 0.36;
DELTA = 0.025;
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
disp('df/dx (A) =')
disp(fx)
num_fx = eval(fx);  A=num_fx;
disp('df/dx (A) =');
disp(num_fx);

disp('df/dxb (B) =')
disp(fxb)
num_fxb = eval(fxb); B=num_fxb;
disp('df/dxb (B) =');
disp(num_fxb);
%%
disp('df/dy (C) =')
disp(fyt)
num_fy = eval(fyt); C=num_fy;
disp('df/dy (C) =');
disp(num_fy);
%%
disp('df/dz (D) =')
disp(fz)
num_fz = eval(fz); D=num_fz;
disp('df/dz (D)=');
disp(num_fz');

disp('df/dx (F) =')
disp(eq1xf)
num_eq1xf = eval(eq1xf); F=num_eq1xf;
disp('df/dx (F) =');
disp(num_eq1xf);

disp('df/dxb (G) =')
disp(eq1x)
num_eq1x = eval(eq1x);  G=num_eq1x;
disp('df/dx (G) =');
disp(num_eq1x);

disp('df/dxb (H) =')
disp(eq1xb)
num_eq1xb = eval(eq1xb); H=num_eq1xb;
disp('df/dx (H) =');
disp(num_eq1xb);

disp('df/dyf (J) =')
disp(eq1yf)
num_eq1yf = eval(eq1yf); J=num_eq1yf;
disp('df/dyf (J) =');
disp(num_eq1yf);

disp('df/dy (K) =')
disp(eq1yt)
num_eq1y = eval(eq1yt);  K=num_eq1y;
disp('df/dy (K) =');
disp(num_eq1y);

disp('df/dzf (L) =')
disp(eq1zf)
num_eq1zf = eval(eq1zf);  L=num_eq1zf;
disp('df/dzf (L) =');
disp(num_eq1zf);

%%
%%  Page 104:  Sec 6.3.1
% make matrix 
% A=[0 -kbar1 0 0]';                         % 4 times1
% B=[0 (1-delta)*kbar1 theta -1]';   % 4 times1
% C=[1 -1 -1/(1-hbar1) 0;...             %  4 times 4
%     ybar -cbar 0 0; ...
%     -1 0 1-theta 0;...
%     1 0 0 -1];
% D=[0 0 1 0]';    % 4 times1
% F=[0];
% G=F;
% H=F;
% J=[0 -1 0 beta*rbar];  % 1 times 4
% K=[0 1 0 0];                    % 1 times 4      
% L=F;
M=[0];
N=[.95];

Cinv=inv(C);

a=F-J*Cinv*A;
b=-(J*Cinv*B-G+K*Cinv*A);
c=-K*Cinv*B+H;

% 2次方程式の解の公式
P1=(-b+sqrt(b^2-4*a*c))/(2*a);
P2=(-b-sqrt(b^2-4*a*c))/(2*a);


if abs(P1)<1
    P=P1;
else
    P=P2;
end
R=-Cinv*(A*P+B);
Q=(J*Cinv*D-L)*N+K*Cinv*D-M;
QD=kron(N',(F-J*Cinv*A))+(J*R+F*P+G-K*Cinv*A);
Q=Q/QD;
S=-Cinv*(A*Q+D);

%% Impulse response

eps = 0.01;
time = 100;

% t=1
   lambda(1)=eps;
   gamma=0.95;
   K(1) = Q*  lambda(1);
   y(:,1)= S* lambda(1);
  
 % t >= 2  
for  t =2:time 
     lambda(t)=  gamma*lambda(t-1);  % 6.10
     K(t)  =  P*K(t-1)+Q*lambda(t);
      y(:,t)=  R* K(t-1)+S*lambda(t);
end

%%
% y =[ Y, C, H, r]

figure('Name','IRF')
l1= plot(lambda,'b-','Linewidth',2);
hold on
   l2=plot(K,'r-.','Linewidth',2);
   l3=plot(y(1,:)','g--','Linewidth',2);
    l4=plot(y(2,:)','m-.','Linewidth',2);
   l5=plot(y(3,:)','c-','Linewidth',2);
   l6=plot(y(4,:)','k:','Linewidth',2);
hold off

xlabel('Time')
set(gca, 'Fontsize',12)
legend([l1, l2, l3, l4,l5,l6],...
    {'lamdda','k','y', 'c' ,'h','r'},'Fontsize',14)
title('Impulse response (log linear)','Fontsize',16)

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
