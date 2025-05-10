function [G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=qzmlre2gensys(A,B,C,c,div)
%% ==============================================================
%% This code transforms the LRE model representation in Cho(2018) into 
%% gensys form and apply gensys algorithm.
%%      Transform A, B and C into g0, g1,c, psi and pi  
%%      where A, B, and C are the coefficients of the LRE model: 
%%          x(t)=A*E[x(t+1)|I(t)]+B x(t-1)+C*z(t),
%% Then this code apply gensys and yields its output.
%% This code use gensys.m, qzdiv and qzswitch.m as nested functions.
%  Input:   A : n by n matrix (forward looking, may well be singular)
%           B : (Optional) n by n matrix (backwardlooking, may well be singular)
%                           Default : B=zeros(n,n);
%           C1 : (Optional) n by m coefficient matrix of the m by 1 vector of
%                exogenous variables, z(t). (Default: C=eye(n))
%           c  : (Optional): This is an optional input of gensys.m
%           div: (Optional): This is an optional input of gensys.m       
%% ==============================================================


%% Transformation
    nn=size(A,1); 
    if nargin==1, B=[]; C=[]; c=[]; div=[]; end
    if nargin==2,       C=[]; c=[]; div=[]; end
    if nargin==3,             c=[]; div=[]; end
    if nargin==4,                   div=[]; end
    
    if isempty(B),      B=zeros(nn,nn); end  
    if isempty(C),      C=eye(nn,nn);   end  
                        mm=size(C,2);
    if isempty(c),      c=zeros(2*nn,1); end

        g0=[eye(nn) -A;eye(nn) zeros(nn,nn)];
        g1=[B zeros(nn,nn); zeros(nn,nn) eye(nn)];
        psi=[C;zeros(nn,mm)];
        pi=[zeros(nn,nn);eye(nn)];
%% Apply gensys
    if isempty(div)
        [G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi);
    else
        [G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi,div);
    end
        
%% The Following gensys.m, qzdiv.m and qzswitch.m are taken from Sims'.    
%% Warning Display is shut off
%% [1] Nested function gensys.m            

function [G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi,div)
% function [G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi,div)
% System given as
%        g0*y(t)=g1*y(t-1)+c+psi*z(t)+pi*eta(t),
% with z an exogenous variable process and eta being endogenously determined
% one-step-ahead expectational errors.  Returned system is
%       y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1) .
% If z(t) is i.i.d., the last term drops out.
% If div is omitted from argument list, a div>1 is calculated.
% eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu(1)=-1 for
% existence only with not-s.c. z; eu=[-2,-2] for coincident zeros.
% By Christopher A. Sims
% Corrected 10/28/96 by CAS
eu=[0;0];
realsmall=1e-6;
fixdiv=(nargin==6);
n=size(g0,1);
[a b q z v]=qz(g0,g1);
if ~fixdiv, div=1.01; end
nunstab=0;
zxz=0;
for i=1:n
% ------------------div calc------------
   if ~fixdiv
      if abs(a(i,i)) > 0
         divhat=abs(b(i,i))/abs(a(i,i));
	 % bug detected by Vasco Curdia and Daria Finocchiaro, 2/25/2004  A root of
	 % exactly 1.01 and no root between 1 and 1.02, led to div being stuck at 1.01
	 % and the 1.01 root being misclassified as stable.  Changing < to <= below fixes this.
         if 1+realsmall<divhat & divhat<=div
            div=.5*(1+divhat);
         end
      end
   end
% ----------------------------------------
   nunstab=nunstab+(abs(b(i,i))>div*abs(a(i,i)));
   if abs(a(i,i))<realsmall & abs(b(i,i))<realsmall
      zxz=1;
   end
end
div ;
nunstab;
if ~zxz
   [a b q z]=qzdiv(div,a,b,q,z);
end
gev=[diag(a) diag(b)];
if zxz
   disp('Coincident zeros.  Indeterminacy and/or nonexistence.')
   eu=[-2;-2];
   % correction added 7/29/2003.  Otherwise the failure to set output
   % arguments leads to an error message and no output (including eu).
   G1=[];C=[];impact=[];fmat=[];fwt=[];ywt=[];gev=[];
   return
end
q1=q(1:n-nunstab,:);
q2=q(n-nunstab+1:n,:);
z1=z(:,1:n-nunstab)';
z2=z(:,n-nunstab+1:n)';
a2=a(n-nunstab+1:n,n-nunstab+1:n);
b2=b(n-nunstab+1:n,n-nunstab+1:n);
etawt=q2*pi;
% zwt=q2*psi;
% branch below is to handle case of no stable roots, which previously (5/9/09)
% quit with an error in that case.
neta = size(pi,2);
if nunstab == 0
  etawt == zeros(0,neta);
  ueta = zeros(0,0);
  deta = zeros(0,0);
  veta = zeros(neta,0);
  bigev = 0;
else
  [ueta,deta,veta]=svd(etawt);
  md=min(size(deta));
  bigev=find(diag(deta(1:md,1:md))>realsmall);
  ueta=ueta(:,bigev);
  veta=veta(:,bigev);
  deta=deta(bigev,bigev);
end
% ------ corrected code, 3/10/04
eu(1) = length(bigev)>=nunstab;
% ------ Code below allowed "existence" in cases where the initial lagged state was free to take on values
% ------ inconsistent with existence, so long as the state could w.p.1 remain consistent with a stable solution
% ------ if its initial lagged value was consistent with a stable solution.  This is a mistake, though perhaps there
% ------ are situations where we would like to know that this "existence for restricted initial state" situation holds.
%% [uz,dz,vz]=svd(zwt);
%% md=min(size(dz));
%% bigev=find(diag(dz(1:md,1:md))>realsmall);
%% uz=uz(:,bigev);
%% vz=vz(:,bigev);
%% dz=dz(bigev,bigev);
%% if isempty(bigev)
%% 	exist=1;
%% else
%% 	exist=norm(uz-ueta*ueta'*uz) < realsmall*n;
%% end
%% if ~isempty(bigev)
%% 	zwtx0=b2\zwt;
%% 	zwtx=zwtx0;
%% 	M=b2\a2;
%% 	for i=2:nunstab
%% 		zwtx=[M*zwtx zwtx0];
%% 	end
%% 	zwtx=b2*zwtx;
%% 	[ux,dx,vx]=svd(zwtx);
%% 	md=min(size(dx));
%% 	bigev=find(diag(dx(1:md,1:md))>realsmall);
%% 	ux=ux(:,bigev);
%% 	vx=vx(:,bigev);
%% 	dx=dx(bigev,bigev);
%% 	existx=norm(ux-ueta*ueta'*ux) < realsmall*n;
%% else
%% 	existx=1;
%% end
% ----------------------------------------------------
% Note that existence and uniqueness are not just matters of comparing
% numbers of roots and numbers of endogenous errors.  These counts are
% reported below because usually they point to the source of the problem.
% ------------------------------------------------------
% branch below to handle case of no stable roots
if nunstab == n
  etawt1 = zeros(0,neta);
  bigev =0;
  ueta1 = zeros(0, 0);
  veta1 = zeros(neta,0);
  deta1 = zeros(0,0);
else
  etawt1 = q1 * pi;
  ndeta1 = min(n-nunstab,neta);
  [ueta1,deta1,veta1]=svd(etawt1);
  md=min(size(deta1));
  bigev=find(diag(deta1(1:md,1:md))>realsmall);
  ueta1=ueta1(:,bigev);
  veta1=veta1(:,bigev);
  deta1=deta1(bigev,bigev);
end
%% if existx | nunstab==0
%%    %disp('solution exists');
%%    eu(1)=1;
%% else
%%     if exist
%%         %disp('solution exists for unforecastable z only');
%%         eu(1)=-1;
%%     %else
%%         %fprintf(1,'No solution.  %d unstable roots. %d endog errors.\n',nunstab,size(ueta1,2));
%%     end
%%     %disp('Generalized eigenvalues')
%%    %disp(gev);
%%    %md=abs(diag(a))>realsmall;
%%    %ev=diag(md.*diag(a)+(1-md).*diag(b))\ev;
%%    %disp(ev)
%% %   return;
%% end
if isempty(veta1)
	unique=1;
else
	loose = veta1-veta*veta'*veta1;
	[ul,dl,vl] = svd(loose);
	nloose = sum(abs(diag(dl)) > realsmall*n);
	unique = (nloose == 0);
end
if unique
   %disp('solution unique');
   eu(2)=1;
else
   %fprintf(1,'Indeterminacy.  %d loose endog errors.\n',nloose);  =======> SHUT DOWN MESSAGE in this code
   %disp('Generalized eigenvalues')
   %disp(gev);
   %md=abs(diag(a))>realsmall;
   %ev=diag(md.*diag(a)+(1-md).*diag(b))\ev;
   %disp(ev)
%   return;
end
tmat = [eye(n-nunstab) -(ueta*(deta\veta')*veta1*deta1*ueta1')'];
G0= [tmat*a; zeros(nunstab,n-nunstab) eye(nunstab)];
G1= [tmat*b; zeros(nunstab,n)];
% ----------------------
% G0 is always non-singular because by construction there are no zeros on
% the diagonal of a(1:n-nunstab,1:n-nunstab), which forms G0's ul corner.
% -----------------------
G0I=inv(G0);
G1=G0I*G1;
usix=n-nunstab+1:n;
C=G0I*[tmat*q*c;(a(usix,usix)-b(usix,usix))\q2*c];
impact=G0I*[tmat*q*psi;zeros(nunstab,size(psi,2))];
fmat=b(usix,usix)\a(usix,usix);
fwt=-b(usix,usix)\q2*psi;
ywt=G0I(:,usix);
% Correction 5/07/2009:  formerly had forgotten to premultiply by G0I
loose = G0I * [etawt1 * (eye(neta) - veta * veta');zeros(nunstab, neta)];
% -------------------- above are output for system in terms of z'y -------
G1=real(z*G1*z');
C=real(z*C);
impact=real(z*impact);
loose = real(z * loose);
% Correction 10/28/96:  formerly line below had real(z*ywt) on rhs, an error.
ywt=z*ywt;

end %% end of gensys.m


%% [2] Nested function qzdiv.m
function [A,B,Q,Z,v] = qzdiv(stake,A,B,Q,Z,v)
%function [A,B,Q,Z,v] = qzdiv(stake,A,B,Q,Z,v)
%
% Takes U.T. matrices A, B, orthonormal matrices Q,Z, rearranges them
% so that all cases of abs(B(i,i)/A(i,i))>stake are in lower right 
% corner, while preserving U.T. and orthonormal properties and Q'AZ' and
% Q'BZ'.  The columns of v are sorted correspondingly.
%
% by Christopher A. Sims
% modified (to add v to input and output) 7/27/00
vin = nargin==6;
if ~vin, v=[]; end;
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
      if vin
         tmp=v(:,k);
         v(:,k)=v(:,k+1);
         v(:,k+1)=tmp;
      end
   end
end        
end %% End of qzdiv.m

function [A,B,Q,Z] = qzswitch(i,A,B,Q,Z)
%function [A,B,Q,Z] = qzswitch(i,A,B,Q,Z)
%
% Takes U.T. matrices A, B, orthonormal matrices Q,Z, interchanges
% diagonal elements i and i+1 of both A and B, while maintaining
% Q'AZ' and Q'BZ' unchanged.  If diagonal elements of A and B
% are zero at matching positions, the returned A will have zeros at both
% positions on the diagonal.  This is natural behavior if this routine is used
% to drive all zeros on the diagonal of A to the lower right, but in this case
% the qz transformation is not unique and it is not possible simply to switch
% the positions of the diagonal elements of both A and B.
 realsmall=sqrt(eps)*10;
%realsmall=1e-3;
a = A(i,i); d = B(i,i); b = A(i,i+1); e = B(i,i+1);
c = A(i+1,i+1); f = B(i+1,i+1);
		% A(i:i+1,i:i+1)=[a b; 0 c];
		% B(i:i+1,i:i+1)=[d e; 0 f];
if (abs(c)<realsmall & abs(f)<realsmall)
	if abs(a)<realsmall
		% l.r. coincident 0's with u.l. of A=0; do nothing
		return
	else
		% l.r. coincident zeros; put 0 in u.l. of a
		wz=[b; -a];
		wz=wz/sqrt(wz'*wz);
		wz=[wz [wz(2)';-wz(1)'] ];
		xy=eye(2);
	end
elseif (abs(a)<realsmall & abs(d)<realsmall)
	if abs(c)<realsmall
		% u.l. coincident zeros with l.r. of A=0; do nothing
		return
	else
		% u.l. coincident zeros; put 0 in l.r. of A
		wz=eye(2);
		xy=[c -b];
		xy=xy/sqrt(xy*xy');
		xy=[[xy(2)' -xy(1)'];xy];
	end
else
	% usual case
	wz = [c*e-f*b, (c*d-f*a)'];
	xy = [(b*d-e*a)', (c*d-f*a)'];
	n = sqrt(wz*wz');
	m = sqrt(xy*xy');
	if m<eps*100
		% all elements of A and B proportional
		return
	end
   wz = n\wz;
   xy = m\xy;
   wz = [wz; -wz(2)', wz(1)'];
   xy = [xy;-xy(2)', xy(1)'];
end
A(i:i+1,:) = xy*A(i:i+1,:);
B(i:i+1,:) = xy*B(i:i+1,:);
A(:,i:i+1) = A(:,i:i+1)*wz;
B(:,i:i+1) = B(:,i:i+1)*wz;
Z(:,i:i+1) = Z(:,i:i+1)*wz;
Q(i:i+1,:) = xy*Q(i:i+1,:);

end %% End of qzswitch.m


end
