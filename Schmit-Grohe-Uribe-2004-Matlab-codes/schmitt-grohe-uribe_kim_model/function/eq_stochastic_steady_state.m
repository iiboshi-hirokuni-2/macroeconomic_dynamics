function [RR] = eq_stochastic_steady_state(x0, hx, hxx, hss, sig )

RR = zeros(size(x0));

nx=size(hx,1);

xs = x0;

for i=1:nx
RR(i,1)= -xs(i,1) + hx(i,:) * xs(:,1) ...
         + 1/2 * xs(:,1)' * squeeze(hxx(i,:,:)) * xs(:,1) + 1/2 * hss(i,1)*sig^2;
end
 x0= xs;

% xf= x0(1:nx,1);
% xs= x0(nx+1:2*nx,1); %zeros(size(x0));
% 
% for i=1:nx
% RR(i,1)= -xf(i,1) + hx(i,:) * xf(:,1) ;
% RR(nx+i,1)= -xs(i,1) + hx(i,:) * xs(:,1) ...
%          + 1/2 * xf(:,1)' * squeeze(hxx(i,:,:)) * xf(:,1) + 1/2 * hss(i,1)*sig^2;
% end
% 
%  x0= [ xf; xs ];



