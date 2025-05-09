

clc

tic

syms THETA DELTA A BETTA GAMMA;
syms k kp h lam lamb;
syms eps;


% Production Function
f = lam*k^THETA*h^(1-THETA);

% budget constraint
c = f + (1-DELTA)*k-kp;

% utility function
% eq1 = log(exp(k)^THETA*exp(h)^(1-THETA)+ (1-DELTA)*exp(k)-exp(kp)) + A*log(1-exp(h));
 eq1 = log(c) + A*log(1-h);

% shock
 eq2 = -lam + GAMMA*lamb + eps ;

% variables
z = [ k,  lam, kp, h];

% Compute the first and second derivatives of f
fz=jacobian(eq1,z);

fzz=jacobian(fz',z);


% set parameters
BETTA = 0.99;
THETA = 0.36;
DELTA = 0.025;
GAMMA = 0.95;
A = 1.72;

%% solve steatdystate
%x=  [kp, y, c h];
x0= [0, 0, 0, 0 ];
fun = @(x)steatdystate(x,para);
x1= fsolve(fun,x0);

k_ss=exp(x1(1)); y_ss=exp(x1(2)); c_ss=exp(x1(3)); h_ss=exp(x1(4));
k=k_ss; kp=k_ss; yy=y_ss; c=c_ss; h=h_ss; 
lam = 1;

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

u=log(yy-DELTA*k)+A*log(1-h);

z = [ k,  lam, kp, h];

 m(1,1)= eval(log(yy-DELTA*k)+A*log(1-h)-z*fz'+z*fzz*z'/2);
 m(1,2:5)= eval(fz/2-z*fzz/2);
 m(2:5,1) = m(1,2:5)';
 m(2:5,2:5)= num_fzz;
 
 disp('M =')
 disp(m);
 
% page 161  
  AA=[1       0    0;
      0       0    0;
      1-GAMMA 0  GAMMA ];

  B=  [0 0;
       1 0;
       0 0 ];

  C = [0; 0;  1 ]; 
 
%% cal 
 R=m(1:3,1:3);
 Q=m(4:5,4:5);
 W=m(1:3,4:5)';

  P= eye(3); 

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
 
 %% Impulse response
 eps = 0.01;
 time = 100;
    yy= zeros(7,time); 
    yy(1,1)=1;
    yy(2,1)=k_ss;
    yy(3,1)=1;
    yy(1:3,1) = yy(1:3,1) + C*eps;
% t >= 2  
 for  t =2:time         
    
       yy(1,t)=1; % constant
       lam= GAMMA*yy(3,t-1)+(1-GAMMA)*1; % lam
      % capital
       k=F(1,:)*yy(1:3,t-1); % k       
       kp=F(1,:)*[1,k,lam]'; % kp
      % labor
      h=F(2,:)*yy(1:3,t-1); % h
      % production function
      y= lam*k^THETA*h^(1-THETA); %y
      % consumption
      c= y + (1-DELTA)*k- kp; % c

      yy(2:7,t)=[k; lam; kp; h; y; c ];
 end
 
 %%
% y =[ 1, k, lam ]
% 
 figure('Name','IRF')
 
 hold on
     l1=plot(log(yy(2,:)'./k_ss),'Linewidth',2);
     l2=plot(yy(3,:)'-lam,'Linewidth',2);
     l3=plot(log(yy(5,:)'./h_ss),'Linewidth',2);
     l4=plot(log(yy(6,:)'./y_ss),'Linewidth',2);
     l5=plot(log(yy(7,:)'./c_ss),'Linewidth',2);
 hold off
% 
 xlabel('Time')
 set(gca, 'Fontsize',12)
 legend([l1, l2, l3, l4, l5],...
     {'k','lambda','h','y','c' },'Fontsize',14)
 title('Impulse response (LQDP)','Fontsize',16)

disp([ 'cal time =' num2str(toc) 'sec' ])



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
