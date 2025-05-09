
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
x0= [0, 0, 0, 0, 0];
fun = @(x)steatdystate(x,para);
x1= fsolve(fun,x0);
x1 = real(x1);

k_ss=exp(x1(1)); y_ss=exp(x1(2)); c_ss=exp(x1(3)); h_ss=exp(x1(4)); r_ss =exp(x1(5));

disp('Steady State :')
disp([ 'k_ss= ', num2str(exp(x1(1)))]);
disp([ 'y_ss= ', num2str(exp(x1(2)))]);
disp([ 'c_ss= ', num2str(exp(x1(3)))]);
disp([ 'h_ss= ', num2str(exp(x1(4)))]);
disp([ 'r_ss= ', num2str(exp(x1(5)))]);
% disp([ 'lam_ss= ', num2str(exp(x1(6)))]);

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

disp([ 'cal time =' num2str(toc) 'sec' ])

%%
function R = steatdystate(x,para)

BETTA=para(1); THETA=para(2); DELTA=para(3); A=para(4);

% [kp, y, c h, r lam];
k=x(1); y=x(2); c=x(3); h=x(4); r =x(5); lam=0;
cp=x(3); rp=x(5); kp=x(1);

h0= 0.583;
B= A*log(1-h0)/h0;

R(1) = exp(h)-1/(1+A/(1-THETA)*(1-DELTA*THETA*BETTA/(1-BETTA*(1-DELTA)))) ;

R(2) = exp(k)-(THETA*BETTA/(1-BETTA*(1-DELTA)))^(1/(1-THETA))*exp(h);

R(3) = exp(y) - DELTA*exp(k) - exp(c);

R(4) = exp(lam)*exp(k)^THETA*exp(h)^(1-THETA) -exp(y) ;

R(5) = THETA*exp(lam)*exp(k)^(THETA-1)*h^(1-THETA) - exp(r);

end

