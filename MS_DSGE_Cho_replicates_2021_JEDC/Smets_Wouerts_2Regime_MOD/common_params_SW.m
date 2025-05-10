%     ctou        ${\delta}$          (long_name='depreciation rate')  
%     clandaw     ${\phi_w}$          (long_name='Gross markup wages')   
%     cg          ${\frac{\bar g}{\bar y}}$     (long_name='steady state exogenous spending share')  
%     curvp       ${\varepsilon_p}$   (long_name='Curvature Kimball aggregator prices')  
% parameters curvw ${\varepsilon_w}$  (long_name='Curvature Kimball aggregator wages')  

% // fixed parameters
ctou=0.03;  %.025;
clandaw=1.5;
cg=0.18;
curvp=10;
curvw=10;

% // estimated parameters initialisation
%     calfa       ${\alpha}$          (long_name='capital share')  
%     csigma      ${\sigma_c}$        (long_name='risk aversion')
%     cfc         ${\phi_p}$          (long_name='fixed cost share')  
%     cgy         ${\rho_{ga}}$       (long_name='Feedback technology on exogenous spending')  

calfa=.24;
cbeta=0.98; %.9995;
csigma=1.5;
cfc=1.5;
cgy=0.51;

%     csadjcost   ${\varphi}$         (long_name='investment adjustment cost')  
%     chabb       ${\lambda}$         (long_name='external habit degree')  
%     cindw       ${\iota_w}$         (long_name='Indexation to past wages')  
%     cprobw      ${\xi_w}$           (long_name='Calvo parameter wages') 
%     csigl       ${\sigma_l}$        (long_name='Frisch elasticity')   
%     cindp       ${\iota_p}$         (long_name='Indexation to past prices')  
%     cprobp      ${\xi_p}$           (long_name='Calvo parameter prices')   
%     czcap       ${\psi}$            (long_name='capacity utilization cost')  
csadjcost= 6.0144;
chabb=    0.6361;    
cprobw=   0.8087;
csigl=    1.9423;
cprobp=   0.6;
cindw=    0.3243;
cindp=    0.47;
czcap=    0.2696;

%     constelab   ${\bar l}$          (long_name='steady state hours')  
%     constepinf  ${\bar \pi}$        (long_name='steady state inflation rate')  
%     constebeta  ${100(\beta^{-1}-1)}$ (long_name='time preference rate in percent')  
%     ctrend      ${\bar \delta}$     (long_name='net growth rate in percent')  
constelab=0;

%% Added by JP to provide full calibration of model before estimation
constepinf=0.7;
constebeta=0.7420;
ctrend=0.3982;

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
%     cmaw        ${\mu_w}$           (long_name='coefficient on MA term wage markup')  
%     cmap        ${\mu_p}$           (long_name='coefficient on MA term price markup')  

% crhoa=    0.9977;
% crhob=    0.5799;
% crhog=    0.9957;
% crhols=   0.9928;
% crhoqs=   0.7165;
% crhoas=1; 
% crhoms=0;
% crhopinf=0;
% crhow=0;
% cmap = 0;
% cmaw  = 0;

%     cinvs       ${d_3}$             (long_name='Unused parameter')  
%     ccs         ${d_4}$             (long_name='Unused parameter')  
%     crdpi       ${r_{\Delta \pi}}$  (long_name='Unused parameter')  
