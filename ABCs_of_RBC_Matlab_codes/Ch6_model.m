
clc

%Define parameters
syms THETA BETTA DELTA A

%Define variables 
syms c  r  y  h  kp lam;  % control variables
syms cf rf yf hf kpf lamf; % future period 
syms cb rb yb hb k  lamb; % state and lagged variables   


%% Set FOC of Model: Page 100 
eq1 = BETTA*(exp(c)/exp(cf))*(exp(rf)+(1-DELTA)) -1 ;

eq2 = (1- THETA)*(1-exp(h))*exp(y)/exp(h)- A*exp(c);

eq3 = exp(y) + (1-DELTA)*exp(k) - exp(kp)- exp(c);

eq4 = exp(lam)*exp(k)^THETA*exp(h)^(1-THETA) -exp(y) ;

eq5 = THETA*exp(y)/exp(k) - exp(r);

%% Create function f
f = [eq1;eq2; eq3; eq4; eq5];

xf = [kpf, yf, cf, hf, rf];
x = [kp, y, c h, r]; % control variables
xb= [k, yb, cb hb, rb]; 
z= lam; % exogenous variables

n = size(f,1); % number of equations
nx = size(x,1); % number of variables

%% Compute the first  derivatives of f
fxf=jacobian(f,xf);
fx=jacobian(f,x);
fxb=jacobian(f,xb);
fz=jacobian(f,z);


% set parameters
BETTA = 0.99;
THETA = 0.36;
DELTA = 0.025;
A = 1.72;
para(1)=BETTA; para(2)=THETA; para(3)=DELTA; para(4)= A;

%% solve steatdystate
%x=  [kp, y, c h, r];
x0= [0, 0, 0, 0, 0,0];
fun = @(x)steatdystate(x,para);
x1= fsolve(fun,x0);

k=x1(1); y=x1(2); c=x1(3); h=x1(4); r =x1(5); lam=x1(6);
cf=x1(3); rf=x1(5); kp=x1(1);

disp('Steady State :')
disp([ 'k_ss= ', num2str(exp(x1(1)))]);
disp([ 'y_ss= ', num2str(exp(x1(2)))]);
disp([ 'c_ss= ', num2str(exp(x1(3)))]);
disp([ 'h_ss= ', num2str(exp(x1(4)))]);
disp([ 'r_ss= ', num2str(exp(x1(5)))]);
disp([ 'lam_ss= ', num2str(exp(x1(6)))]);
disp(' '); 

%% 

disp('df/dxf (F) =')
disp(fxf)
num_fxf = eval(fxf); 
disp('df/dxf (F) =');
disp(num_fxf);

disp('df/dx (G) =')
disp(fx)
num_fx = eval(fx); 
disp('df/dx (G) =');
disp(num_fx);

disp('df/dxb (H) =')
disp(fxb)
num_fxb = eval(fxb); 
disp('df/dxb (H) =');
disp(num_fxb);

disp('df/dz (M) =')
disp(fz)
num_fz = eval(fz); 
disp('df/dz (M)=');
disp(num_fz');

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

