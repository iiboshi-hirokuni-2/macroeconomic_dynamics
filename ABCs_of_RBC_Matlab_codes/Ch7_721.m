
clc

syms THETA DELTA A BETTA;
syms k kp h;


% Production Function
f = k^THETA*h^(1-THETA);

% budget constraint
c = f + (1-DELTA)*k-kp;

% utility function
% eq1 = log(exp(k)^THETA*exp(h)^(1-THETA)+ (1-DELTA)*exp(k)-exp(kp)) + A*log(1-exp(h));
 eq1 = log(c) + A*log(1-h);

% variables
z = [ k,  kp, h];

% Compute the first and second derivatives of f
% fz=jacobian(eq1,z);
fz=jacobian(eq1,z);

fzz=jacobian(fz',z);


% set parameters
BETTA = 0.99;
THETA = 0.36;
DELTA = 0.025;
A = 1.72;

%% solve steatdystate
%x=  [kp, y, c h];
x0= [0, 0, 0, 0 ];
fun = @(x)steatdystate(x,para);
x1= fsolve(fun,x0);

k_ss=exp(x1(1)); y_ss=exp(x1(2)); c_ss=exp(x1(3)); h_ss=exp(x1(4));
k=k_ss; kp=k_ss; y=y_ss; c=c_ss; h=h_ss; 

disp('Steady State :')
disp([ 'k_ss= ', num2str(exp(x1(1)))]);
disp([ 'y_ss= ', num2str(exp(x1(2)))]);
disp([ 'c_ss= ', num2str(exp(x1(3)))]);
disp([ 'h_ss= ', num2str(exp(x1(4)))]);
disp(' '); 

%%

disp(' ');
disp('df/dz  =')
disp(fz')
num_fz = eval(fz); 
disp('df/dz =');
disp(num_fz);

disp(' ');
disp('df/dz^2 =')
disp(fzz)
num_fzz = eval(fzz)/2; B=num_fzz;
disp('df/dzz^2 =');
disp(num_fzz);

u=log(y-DELTA*k)+A*log(1-h);

m(1,1)= eval(log(y-DELTA*k)+A*log(1-h)-[k,kp,h]*fz'+[k,kp,h]*fzz*[k,kp,h]'/2);
m(1,2:4)= eval(fz/2-[k,kp,h]*fzz/2);
m(2:4,1) = m(1,2:4)';
m(2:4,2:4)= num_fzz;

disp('M =')
disp(m);

AA=[1 0
    0 0];
B=[0 0
    1 0];


%% cal 
R=m(1:2,1:2);
Q=m(3:4,3:4);
W=m(1:2,3:4)';
P=[1 0
    0 1];
for i=1:1000
    zinv=inv(Q+BETTA*B'*P*B);
    z2=BETTA*AA'*P*B+W';
    P=R+BETTA*AA'*P*AA-z2*zinv*z2';
end

F=-zinv*(W+BETTA*B'*P*AA);

disp('P=')
disp(P)
disp('F=')
disp(F)

%%
function R = steatdystate(x,para)

BETTA=para(1); THETA=para(2); DELTA=para(3); A=para(4);

% [kp, y, c h, r lam];
k=x(1); y=x(2); c=x(3); h=x(4); 
cp=x(3);  kp=x(1);

R(1) = exp(h)-1/(1+A/(1-THETA)*(1-DELTA*THETA*BETTA/(1-BETTA*(1-DELTA)))) ;

R(2) = exp(k)-(THETA*BETTA/(1-BETTA*(1-DELTA)))^(1/(1-THETA))*exp(h);

R(3) = exp(y) - exp(k)^THETA*exp(h)^(1-THETA) ;

R(4) = exp(y) - exp(c) + (1-DELTA)*exp(k) - exp(k);


end
