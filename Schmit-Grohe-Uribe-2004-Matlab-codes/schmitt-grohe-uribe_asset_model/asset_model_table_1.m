


[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,...
    fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,...
    fxy,fxxp,fxx,eta]=asset_model(P); 

%First-order approximation of Policy function
[gx,hx] = gx_hx(fy,fx,fyp,fxp);

%Second-order approximation of Policy function
[gxx,hxx] = gxx_hxx(fx,fxp,fy,fyp,fypyp,fypy,fypxp,...
                   fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,...
                   fxpxp,fxpx,fxyp,fxy,fxxp,fxx,hx,gx); 

[gss,hss] = gss_hss(fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,...
                    fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,...
                    fxy,fxxp,fxx,hx,gx,gxx,eta);
load('p_bar.mat')

disp('Table 1: 2nd Order approx.')
disp(['THETA = ',  num2str(P.THETA) ]);
disp(['RHO = ',  num2str(P.RHO) ]);
disp(['sig = ',  num2str(P.SIG) ]);
disp(['f0 = ',  num2str(p+1/2*gss*sig) ]);
disp(['f1 = ',  num2str(gx) ]);
disp(['f2 = ',  num2str(gxx) ]);

% disp(' ')
% disp('Table 1: Fixed-point ')
% gx_fixed = P.THETA*P.RHO*BETTA*exp(P.THETA*a)...
%         /(1-BETTA*exp(P.THETA*a))/(1-BETTA*exp(P.THETA*a)*P.RHO) ;
% disp(['f1 = ',  num2str(gx_fixed) ]);




