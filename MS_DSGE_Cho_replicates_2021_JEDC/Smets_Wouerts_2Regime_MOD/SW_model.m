
% disp( '  '  )
% disp('============================================')
% disp('           start  Smets & Wouters (AER 2007) DSGE model ')
% disp('============================================')

% var labobs      ${lHOURS}$      (long_name='log hours worked') 
%     robs        ${FEDFUNDS}$    (long_name='Federal funds rate') 
%     pinfobs     ${dlP}$         (long_name='Inflation') 
%     dy          ${dlGDP}$       (long_name='Output growth rate') 
%     dc          ${dlCONS}$      (long_name='Consumption growth rate') 
%     dinve       ${dlINV}$       (long_name='Investment growth rate') 
%     dw          ${dlWAG}$       (long_name='Wage growth rate') 
%     ewma        ${\eta^{w,aux}}$ (long_name='Auxiliary wage markup moving average variable')  
%     epinfma     ${\eta^{p,aux}}$ (long_name='Auxiliary price markup moving average variable')  
%     zcapf       ${z^{flex}}$    (long_name='Capital utilization rate flex price economy') 
%     rkf         ${r^{k,flex}}$  (long_name='rental rate of capital flex price economy') 
%     kf          ${k^{s,flex}}$  (long_name='Capital services flex price economy') 
%     pkf         ${q^{flex}}$    (long_name='real value of existing capital stock flex price economy')  
%     cf          ${c^{flex}}$    (long_name='Consumption flex price economy') 
%     invef       ${i^{flex}}$    (long_name='Investment flex price economy') 
%     yf          ${y^{flex}}$    (long_name='Output flex price economy') 
%     labf        ${l^{flex}}$    (long_name='hours worked flex price economy') 
%     wf          ${w^{flex}}$    (long_name='real wage flex price economy') 
%     rrf         ${r^{flex}}$    (long_name='real interest rate flex price economy')
%     mc          ${\mu_p}$       (long_name='gross price markup') 
%     zcap        ${z}$           (long_name='Capital utilization rate') 
%     rk          ${r^{k}}$       (long_name='rental rate of capital') 
%     k           ${k^{s}}$       (long_name='Capital services') 
%     pk          ${q}$           (long_name='real value of existing capital stock') 
%     c           ${c}$           (long_name='Consumption')
%     inve        ${i}$           (long_name='Investment')
%     y           ${y}$           (long_name='Output')
%     lab         ${l}$           (long_name='hours worked')
%     pinf        ${\pi}$         (long_name='Inflation')
%     w           ${w}$           (long_name='real wage')
%     r           ${r}$           (long_name='nominal interest rate')
%     a           ${\varepsilon_a}$       (long_name='productivity process')
%     b           ${c_2*\varepsilon_t^b}$ (long_name='Scaled risk premium shock')
%     g           ${\varepsilon^g}$       (long_name='Exogenous spending')
%     qs          ${\varepsilon^i}$       (long_name='Investment-specific technology')
%     ms          ${\varepsilon^r}$       (long_name='Monetary policy shock process') 
%     spinf       ${\varepsilon^p}$       (long_name='Price markup shock process')
%     sw          ${\varepsilon^w}$       (long_name='Wage markup shock process')
%     kpf         ${k^{flex}}$            (long_name='Capital stock flex price economy') 
%     kp          ${k}$           (long_name='Capital stock') 
%     ;    
%  
% varexo ea       ${\eta^a}$      (long_name='productivity shock')
%     eb          ${\eta^b}$      (long_name='Investment-specific technology shock')
%     eg          ${\eta^g}$      (long_name='Spending shock')
%     eqs         ${\eta^i}$      (long_name='Investment-specific technology shock')
%     em          ${\eta^m}$      (long_name='Monetary policy shock')
%     epinf       ${\eta^{p}}$    (long_name='Price markup shock')  
%     ew          ${\eta^{w}}$    (long_name='Wage markup shock')  
%         ;  
%  
% parameters curvw ${\varepsilon_w}$  (long_name='Curvature Kimball aggregator wages')  
%     cgy         ${\rho_{ga}}$       (long_name='Feedback technology on exogenous spending')  
%     curvp       ${\varepsilon_p}$   (long_name='Curvature Kimball aggregator prices')  
%     constelab   ${\bar l}$          (long_name='steady state hours')  
%     constepinf  ${\bar \pi}$        (long_name='steady state inflation rate')  
%     constebeta  ${100(\beta^{-1}-1)}$ (long_name='time preference rate in percent')  
%     cmaw        ${\mu_w}$           (long_name='coefficient on MA term wage markup')  
%     cmap        ${\mu_p}$           (long_name='coefficient on MA term price markup')  
%     calfa       ${\alpha}$          (long_name='capital share')  
%     czcap       ${\psi}$            (long_name='capacity utilization cost')  
%     csadjcost   ${\varphi}$         (long_name='investment adjustment cost')  
%     ctou        ${\delta}$          (long_name='depreciation rate')  
%     csigma      ${\sigma_c}$        (long_name='risk aversion')  
%     chabb       ${\lambda}$         (long_name='external habit degree')  
%     ccs         ${d_4}$             (long_name='Unused parameter')  
%     cinvs       ${d_3}$             (long_name='Unused parameter')  
%     cfc         ${\phi_p}$          (long_name='fixed cost share')  
%     cindw       ${\iota_w}$         (long_name='Indexation to past wages')  
%     cprobw      ${\xi_w}$           (long_name='Calvo parameter wages')   
%     cindp       ${\iota_p}$         (long_name='Indexation to past prices')  
%     cprobp      ${\xi_p}$           (long_name='Calvo parameter prices')   
%     csigl       ${\sigma_l}$        (long_name='Frisch elasticity')   
%     clandaw     ${\phi_w}$          (long_name='Gross markup wages')   
%     crdpi       ${r_{\Delta \pi}}$  (long_name='Unused parameter')  
%     crpi        ${r_{\pi}}$         (long_name='Taylor rule inflation feedback') 
%     crdy        ${r_{\Delta y}}$    (long_name='Taylor rule output growth feedback') 
%     cry         ${r_{y}}$           (long_name='Taylor rule output level feedback') 
%     crr         ${\rho}$            (long_name='interest rate persistence')  
%     crhoa       ${\rho_a}$          (long_name='persistence productivity shock')  
%     crhoas      ${d_2}$             (long_name='Unused parameter')  
%     crhob       ${\rho_b}$          (long_name='persistence risk premium shock')  
%     crhog       ${\rho_g}$          (long_name='persistence spending shock')  
%     crhols      ${d_1}$             (long_name='Unused parameter')  
%     crhoqs      ${\rho_i}$          (long_name='persistence risk premium shock')  
%     crhoms      ${\rho_r}$          (long_name='persistence monetary policy shock')  
%     crhopinf    ${\rho_p}$          (long_name='persistence price markup shock')  
%     crhow       ${\rho_w}$          (long_name='persistence wage markup shock')  
%     ctrend      ${\bar \delta}$     (long_name='net growth rate in percent')  
%     cg          ${\frac{\bar g}{\bar y}}$     (long_name='steady state exogenous spending share')  
%     ;

% // flexible economy
syms yf yf1;

% endogenous variables
syms  mc zcap rk k pk c inve y lab pinf w r kp;   %a b g qs ms   spinf sw kpf kp;   
        %ewma epinfma zcapf rkf kf pkf cf  invef yf labf wf rrf
% lag variables
syms  mc1 zcap1 rk1 k pk1 c1 inve1 y1 lab1 pinf1 w1 r1 kp1 ;
% forward variables
syms  mc0 zcap0 rk0 k0 pk0 c0 inve0 y0 lab0 pinf0 w0 r0 kp0;

% exogenous variables
syms a b g qs ms   spinf sw kpf ;  % ea eb eg eqs em epinf ew;  

% parameters
syms  curvw cgy curvp constelab constepinf constebeta cmaw cmap calfa ...
    czcap csadjcost ctou csigma chabb ccs cinvs cfc ...
    cindw cprobw cindp cprobp csigl clandaw  ...
    crdpi crpi crdy cry crr  ...
    crhoa crhoas crhob crhog crhols crhoqs crhoms crhopinf crhow  ...
    ctrend cg ;

%% make MATRIX
Y_f =  transpose([ mc0 zcap0 rk0 k0 pk0 c0 inve0 y0 lab0 pinf0 w0 r0 kp0  ]);
Y_t = transpose([ mc zcap rk k pk c inve y lab pinf w r kp ]);  
Y_b = transpose([ mc1 zcap1 rk1 k pk1 c1 inve1 y1 lab1 pinf1 w1 r1  kp1  ]);
Epsilon_t = transpose([  a b g qs ms   spinf sw ]);

cpie=1+constepinf/100;
cgamma=1+ctrend/100 ;
cbeta=1/(1+constebeta/100);

clandap=cfc;
cbetabar=cbeta*cgamma^(-csigma);
cr=cpie/(cbeta*cgamma^(-csigma));
crk=(cbeta^(-1))*(cgamma^csigma) - (1-ctou);
cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
%//cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*((cbeta^(-1))*(cgamma^csigma) - (1-ctou))^calfa))^(1/(1-calfa));
cikbar=(1-(1-ctou)/cgamma);
cik=(1-(1-ctou)/cgamma)*cgamma;
clk=((1-calfa)/calfa)*(crk/cw);
cky=cfc*(clk)^(calfa-1);
ciy=cik*cky;
ccy=1-cg-cik*cky;
crkky=crk*cky;
cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
cwly=1-crk*cky;

conster=(cr-1)*100;

%% Model
%% // sticky price - wage economy
%  eq1        [name='FOC labor with mpl expressed as function of rk and w, SW Equation (9)']
% mc =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;
	eq1= - mc +  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;

%	eq2      [name='FOC capacity utilization, SW Equation (7)']
%     zcap =  (1/(czcap/(1-czcap)))* rk ; 
    eq2 =- zcap+  (1/(czcap/(1-czcap)))* rk ;

 %  eq3        [name='Firm FOC capital, SW Equation (11)']
  %     rk =  w+lab-k ;
	 eq3 =  -rk + w+lab-k ;

 %  eq4       [name='Definition capital services, SW Equation (6)']
 %   k =  kp(-1)+zcap ;
	 eq4 =  -k + kp1+zcap ;

  %   eq5     [name='Investment Euler Equation, SW Equation (3)']
  %  inve = (1/(1+cbetabar*cgamma))* (inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ;
	eq5 =  - inve + (1/(1+cbetabar*cgamma))* (inve1 + cbetabar*cgamma*inve0+(1/(cgamma^2*csadjcost))*pk ) +qs ;

%     eq6         [name='Arbitrage equation value of capital, SW Equation (4)']
%   pk = -r+pinf(1)-0*b 
%                 + (1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b 
%                 + (crk/(crk+(1-ctou)))*rk(1) 
%                 + ((1-ctou)/(crk+(1-ctou)))*pk(1) ;
    eq6 =  - pk  -r+pinf0-0*b ...
                + (1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b ...
                + (crk/(crk+(1-ctou)))*rk(1) ...
                + ((1-ctou)/(crk+(1-ctou)))*pk0 ;

%   eq7	      [name='Consumption Euler Equation, SW Equation (2)']
%   c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) 
%                 + (1/(1+chabb/cgamma))*c(+1) 
%                 +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) 
%                 - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) +b ;
	eq7 =  -c + (chabb/cgamma)/(1+chabb/cgamma)*c1 ...
                + (1/(1+chabb/cgamma))*c0 ...
                +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab0 ) ...
                - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf0 + 0*b) +b ;

%   eq8	      [name='Aggregate Resource Constraint, SW Equation (1)']
%          y = ccy*c+ciy*inve+g  +  1*crkky*zcap ;
     eq8 =  -y + ccy*c+ciy*inve+g  +  1*crkky*zcap ;

 %   eq9         [name='Aggregate Production Function, SW Equation (5)']
 %        y = cfc*( calfa*k+(1-calfa)*lab +a );
    eq9 =  -y + cfc*( calfa*k+(1-calfa)*lab +a );

 %   eq10         [name='New Keynesian Phillips Curve, SW Equation (10)']
 %     pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1) 
 %                +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf ; 
  eq10 =  -  pinf + (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf0 +cindp*pinf1 ...
               +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf ; 

%   eq11	      [name='Wage Phillips Curve, SW Equation (13), with (12) plugged for mu_w']
%  w =  (1/(1+cbetabar*cgamma))*w(-1)
%                +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)
%                +(cindw/(1+cbetabar*cgamma))*pinf(-1)
%                -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf
%                +(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)
%                +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*
%                (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w) 
%                + 1*sw ;
  eq11=  - w + (1/(1+cbetabar*cgamma))*w1...
               +(cbetabar*cgamma/(1+cbetabar*cgamma))*w0...
               +(cindw/(1+cbetabar*cgamma))*pinf1...
               -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf...
               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf0...
               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*...
               (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c1 -w) ...
               + 1*sw ;

   %   eq12     [name='Taylor rule, SW Equation (14)']
   %    r =  crpi*(1-crr)*pinf
   %                +cry*(1-crr)*(y-yf)     
   %                +crdy*(y-yf-y(-1)+yf(-1))
   %                +crr*r(-1)
   %                +ms  ;
	 eq12 =  - r +  crpi*(1-crr)*pinf...
               +cry*(1-crr)*(y-yf)     ...
               +crdy*(y-yf-y1+yf1)...
               +crr*r1 ...
               +ms  ;           

   %   eq13       [name='Law of motion for capital, SW Equation (8) (see header notes)']    
  %      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;
	 eq13 =  - kp + (1-cikbar)*kp1+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;      
     
    %% strucural shocks       
%   %   eq        [name='Law of motion for productivity']              
% 	      a = crhoa*a(-1)  + ea;
%    %   eq       [name='Law of motion for risk premium']              
% 	      b = crhob*b(-1) + eb;
%   %   eq        [name='Law of motion for spending process']              
% 	      g = crhog*(g(-1)) + eg + cgy*ea;
% %   eq	      [name='Law of motion for investment specific technology shock process']              
%           qs = crhoqs*qs(-1) + eqs;
%   %   eq        [name='Law of motion for monetary policy shock process']              
% 	      ms = crhoms*ms(-1) + em;
%  %   eq         [name='Law of motion for price markup shock process']              
% 	      spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1);
% 	      epinfma=epinf;
%    %   eq       [name='Law of motion for wage markup shock process']              
% 	      sw = crhow*sw(-1) + ewma - cmaw*ewma(-1) ;
% 	          ewma=ew; 
   

%% Set Sims Forms
System_of_Eq = [  ...
                eq1;  eq2;  eq3;  eq4;  eq5;...
                eq6;  eq7;  eq8;  eq9;  eq10;...
                eq11; eq12; eq13; ...
                ];

neq  = length(System_of_Eq);
% CC = zeros(neq,1);

% chek if all f evaluated at ss value are zero
% check = eval(System_of_Eq)

