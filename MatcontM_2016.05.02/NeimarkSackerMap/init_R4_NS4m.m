function [x0,v0]= init_R4_NS4m(mapfile,eps,x,p,ap,n,varargin)
%
% Initializes a Neimark-Sacker continuation from a R4 point if possible.
%
global cds nsmds
% check input
if size(ap,2)~=2
    errordlg('Two active parameters are needed for a Neimark-Sacker bifurcation curve continuation');
end
v0=[];
% initialize nsmds
nsmds.mapfile = mapfile;
func_handles = feval(nsmds.mapfile);
nsmds.func = func_handles{2};
nsmds.Jacobian  = func_handles{3};
nsmds.JacobianP = func_handles{4};
nsmds.Hessians  = func_handles{5};
nsmds.HessiansP = func_handles{6};
nsmds.Der3      = func_handles{7};
nsmds.Niterations=4*n;
siz = size(func_handles,2);
if siz > 9
    j=1;
    for i=10:siz
        nsmds.user{j}= func_handles{i};
        j=j+1;
    end
end
nsmds.nphase = size(x,1);
nsmds.ActiveParams = ap;
nsmds.P0 = p;
if size(varargin,1)>0,nsmds.BranchParams=varargin{1};else nsmds.BranchParams=[];end 
cds.curve = @neimarksackermap;
cds.ndim = length(x)+2; 
curvehandles = feval(cds.curve);
cds.curve_func = curvehandles{1};
cds.curve_options = curvehandles{3};
%-----Defining Symbolic derivatives-----
  symjac  = ~isempty(nsmds.Jacobian);
  symhes  = ~isempty(nsmds.Hessians);
  symDer3 = ~isempty(nsmds.Der3);
  symord = 0; 
  if symjac, symord = 1; end
  if symhes, symord = 2; end
  if symDer3, symord = 3; end
  cds.options.SymDerivative = symord;
  symjacp  = ~isempty(nsmds.JacobianP); 
  symhessp = ~isempty(nsmds.HessiansP); 
  symordp = 0;
  if symjacp,  symordp = 1; end
  if symhessp, symordp = 2; end
  cds.options.SymDerivativeP = symordp;
%---Branch Switching Algorithm----
    p = n2c(p);nphase = size(x,1);
    A = nsmjac(x,p,n);					
    [X,D] = eig(A-exp(sqrt(-1)*pi/2)*eye(nphase));
    [Y,i] = min(abs(diag(D)));
    vext = X(:,i)/norm(X(:,i));
    [X,D] = eig(A'+exp(sqrt(-1)*pi/2)*eye(nphase));
    [Y,i] = min(abs(diag(D)));
    wext = X(:,i)/(vext'*X(:,i));
    [A0,D0] = nf_R4m(nsmds.func,nsmds.Jacobian,nsmds.Hessians,nsmds.Der3,A,vext,wext,nphase,x,p,n);
    if (1<A0(1)^2) || (abs(A0(2))< (1+A0(1)^2)/sqrt(1-A0(1)^2))
      printconsole('Switching not possible!\n');
      return;
    end 
    hessIncrement =(cds.options.Increment)^(3.0/4.0);
    global T1global T2global
    if (cds.options.SymDerivative >= 2)
      T1global=tens1(nsmds.func,nsmds.Jacobian,x,p,n);
      T2global=tens2(nsmds.func,nsmds.Hessians,x,p,n);
    end
    A1 = nsmjacp(x,p,n);   							%jacobianp
    temp = (eye(nphase)-A)\A1;							%temp=(I-A)^{INV}*J1
    s1=[1;0];s2=[0;1];								%define standard vectors
    xit=zeros(nphase,n);xit(:,1)=x;
    AA=zeros(nphase,nphase,n);
    AA(:,:,1)=nsmjac(x,p,1);
    xx1=x;
    for m=2:n
     xx1=feval(nsmds.func,0,xx1,p{:});
     xit(:,m)=xx1;
     AA(:,:,m)=nsmjac(xx1,p,1);
    end
    test1 = nshesspvect(xit,p,vext,AA,n)*s1;  					% A1(q,s1) 
    test1 = test1 + multilinear2(nsmds.func,vext,temp*s1,x,p,n,hessIncrement); 	% +B(q,temp*s1))
    gamma1= wext'*test1;
    test1 = nshesspvect(xit,p,vext,AA,n)*s2;					% A1(q,s2) 
    test1 = test1 + multilinear2(nsmds.func,vext,temp*s2,x,p,n,hessIncrement); 	% +B(q,temp*s2)
    gamma2= wext'*test1;
    vv = conj([gamma2;-gamma1])/(gamma1*conj(gamma2)-gamma2*conj(gamma1));
    VV = 2*[ real(vv) -imag(vv) ];
    dir= VV*inv([0 4;-4 0])*[-2*A0(1) ; -(A0(2)-sign(A0(2))*sqrt(1-A0(1)*A0(1)))]; % parameter direction
    phi0 = sign(A0(2))*acos(A0(1))/4.0;
    gamma = exp(sqrt(-1.0)*angle(D0)/4)/sqrt(abs(D0));
    xx = 2*real(exp(sqrt(-1.0)*phi0)*gamma*vext);					   % phase direction
    x0 = [x + sqrt(eps)*xx ;nsmds.P0(ap) + eps*dir];			   % predicted point
    clear T1global T2global
%-----End of branch prediction-----------------
[x1,p] = rearr(x0); p = n2c(p);
cds.curve_jacobian =curvehandles{4};
cds.curve_hessians = curvehandles{5};
cds.options = feval(cds.curve_options); 
cds.options = contset(cds.options,'Increment',1e-5);
n = nsmds.Niterations;
nphase = size(x1,1);
nap = length(nsmds.ActiveParams);
% calculate eigenvalues and eigenvectors
jac = nsmjac(x1,p,n);
[V,D] = eig(jac);
% find pair of complex eigenvalues
d = diag(D);
idx1=0;idx2=0;
    for s=1:nphase
       for j=s+1:nphase
            if (abs(1-d(s)*d(j))<0.01)
               idx1=s; 
               idx2=j;
           end
            end
        end
        
if idx1==0
       printconsole('Neutral saddle\n');
        result='Neutral saddle'
       return;
end 
V=V(:,idx1);
[Q,R,E] = qr([real(V(:,idx1)) imag(V(:,idx1))]);
nsmds.borders.v = Q(:,1:2);
[V,D] = eig(jac');
% find pair of complex eigenvalues
d = diag(D);
i=1;
for s=1:nphase
       for j=s+1:nphase
               if (abs(1-d(s)*d(j))<0.01)
               idx1=s;
               idx2=j;
           end
            end
        end
        
if idx1==0
       printconsole('Neutral saddle\n');
        result='Neutral saddle'
       return;
end 
temp=idx1; 
if d(idx1)<0
    idx1=idx2;
    idx2=temp;
end    
[Q,R,E] = qr([real(V(:,idx1)) imag(V(:,idx1))]);
nsmds.borders.w = Q(:,1:2);
k  = real(d(idx1));
x0 = [x0;k];
% calculate eigenvalues
% ERROR OR WARNING
RED  = jac*jac-2*k*jac+eye(nsmds.nphase); 
jacp = nsmjacp(x1,p,n);
A = [jac-eye(nsmds.nphase)  jacp zeros(nsmds.nphase,1)]; 
[Q,R] = qr(A');
Bord  = [RED nsmds.borders.w;nsmds.borders.v' zeros(2)];
bunit = [zeros(nsmds.nphase,2);eye(2)];
vext  = Bord\bunit;
wext  = Bord'\bunit;


AA=zeros(nphase,nphase,n);
xit=zeros(nphase,n);xit(:,1)=x1;
AA(:,:,1)=nsmjac(x1,p,1);
xx1=x1;
for m=2:n
   xx1=feval(nsmds.func,0,xx1,p{:});
   xit(:,m)=xx1;
   AA(:,:,m)=nsmjac(xx1,p,1);
end


wext1=wext(1:nsmds.nphase,1)'*jac;
vext1=vext(1:nsmds.nphase,1);
gx1=nsvecthessvect(xit,p,vext1,wext1,AA,n);

wext2=wext(1:nsmds.nphase,1)';
vext2=jac*vext(1:nsmds.nphase,1);
gx2=nsvecthessvect(xit,p,vext2,wext2,AA,n);

wext3=-2.0*k*wext(1:nsmds.nphase,1)';
vext3=vext(1:nsmds.nphase,1);
gx3=nsvecthessvect(xit,p,vext3,wext3,AA,n);
gxx1=gx1+gx2+gx3;

wext12=wext(1:nsmds.nphase,1)'*jac;
vext12=vext(1:nsmds.nphase,2);
gx12=nsvecthessvect(xit,p,vext12,wext12,AA,n);

wext22=wext(1:nsmds.nphase,1)';
vext22=jac*vext(1:nsmds.nphase,2);
gx22=nsvecthessvect(xit,p,vext22,wext22,AA,n);

wext32=-2.0*k*wext(1:nsmds.nphase,1)';
vext32=vext(1:nsmds.nphase,2);
gx32=nsvecthessvect(xit,p,vext32,wext32,AA,n);
gxx2=gx12+gx22+gx32;
%
wext31=wext(1:nsmds.nphase,2)'*jac;
vext31=vext(1:nsmds.nphase,1);
gx31=nsvecthessvect(xit,p,vext31,wext31,AA,n);

wext32=wext(1:nsmds.nphase,2)';
vext32=jac*vext(1:nsmds.nphase,1);
gx32=nsvecthessvect(xit,p,vext32,wext32,AA,n);

wext33=-2.0*k*wext(1:nsmds.nphase,2)';
vext33=vext(1:nsmds.nphase,1);
gx33=nsvecthessvect(xit,p,vext33,wext33,AA,n);
gxx3=gx31+gx32+gx33;
%
wext41=wext(1:nsmds.nphase,2)'*jac;
vext41=vext(1:nsmds.nphase,2);
gx41=nsvecthessvect(xit,p,vext41,wext41,AA,n);

wext42=wext(1:nsmds.nphase,2)';
vext42=jac*vext(1:nsmds.nphase,2);
gx42=nsvecthessvect(xit,p,vext42,wext42,AA,n);

wext43=-2.0*k*wext(1:nsmds.nphase,2)';
vext43=vext(1:nsmds.nphase,2);
gx43=nsvecthessvect(xit,p,vext43,wext43,AA,n);
gxx4=gx41+gx42+gx43;
%
for i = 1:nsmds.nphase
    gx(1,i)=gxx1(:,i);
    gx(2,i)=gxx2(:,i);
    gx(3,i)=gxx3(:,i);
    gx(4,i)=gxx4(:,i);
    
end

gk(1,1) =2*wext(1:nsmds.nphase,1)'*jac*vext(1:nsmds.nphase,1);
gk(2,1) =2*wext(1:nsmds.nphase,1)'*jac*vext(1:nsmds.nphase,2);
gk(3,1) =2*wext(1:nsmds.nphase,2)'*jac*vext(1:nsmds.nphase,1);
gk(4,1) =2*wext(1:nsmds.nphase,2)'*jac*vext(1:nsmds.nphase,2);

wext1=wext(1:nsmds.nphase,1)'*jac;
vext1=vext(1:nsmds.nphase,1);
gx1=nsvecthesspvect(xit,p,vext1,wext1,AA,n);

wext2=wext(1:nsmds.nphase,1)';
vext2=jac*vext(1:nsmds.nphase,1);
gx2=nsvecthesspvect(xit,p,vext2,wext2,AA,n);

wext3=-2.0*k*wext(1:nsmds.nphase,1)';
vext3=vext(1:nsmds.nphase,1);
gx3=nsvecthesspvect(xit,p,vext3,wext3,AA,n);
gp1=gx1+gx2+gx3;

wext12=wext(1:nsmds.nphase,1)'*jac;
vext12=vext(1:nsmds.nphase,2);
gx12=nsvecthesspvect(xit,p,vext12,wext12,AA,n);

wext22=wext(1:nsmds.nphase,1)';
vext22=jac*vext(1:nsmds.nphase,2);
gx22=nsvecthesspvect(xit,p,vext22,wext22,AA,n);

wext33=-2.0*k*wext(1:nsmds.nphase,1)';
vext33=vext(1:nsmds.nphase,2);
gx32=nsvecthesspvect(xit,p,vext33,wext33,AA,n);
gp2=gx12+gx22+gx32;

wext31=wext(1:nsmds.nphase,2)'*jac;
vext31=vext(1:nsmds.nphase,1);
gx31=nsvecthesspvect(xit,p,vext31,wext31,AA,n);

wext32=wext(1:nsmds.nphase,2)';
vext32=jac*vext(1:nsmds.nphase,1);
gx32=nsvecthesspvect(xit,p,vext32,wext32,AA,n);

wext33=-2.0*k*wext(1:nsmds.nphase,2)';
vext33=vext(1:nsmds.nphase,1);
gx33=nsvecthesspvect(xit,p,vext33,wext33,AA,n);
gp3=gx31+gx32+gx33;

wext41=wext(1:nsmds.nphase,2)'*jac;
vext41=vext(1:nsmds.nphase,2);
gx41=nsvecthesspvect(xit,p,vext41,wext41,AA,n);

wext42=wext(1:nsmds.nphase,2)';
vext42=jac*vext(1:nsmds.nphase,2);
gx42=nsvecthesspvect(xit,p,vext42,wext42,AA,n);

wext43=-2.0*k*wext(1:nsmds.nphase,2)';
vext43=vext(1:nsmds.nphase,2);
gx43=nsvecthesspvect(xit,p,vext43,wext43,AA,n);
gp4=gx41+gx42+gx43;

for i = 1:nap
    gp(1,i)=gp1(:,i);
    gp(2,i)=gp2(:,i);
    gp(3,i)=gp3(:,i);
    gp(4,i)=gp4(:,i);
    
end
A = [A;gx gp gk]*Q;
Jres = A(1+nsmds.nphase:end,1+nsmds.nphase:end)';
[Q,R,E] = qr(Jres');
index = [1 1;1 2;2 1;2 2];
[I,J] = find(E(:,1:2));
nsmds.index1 = index(I(1),:);
nsmds.index2 = index(I(2),:);
rmfield(cds,'options');

% ---------------------------------------------------------------
function [x,p] = rearr(x0)
% [x,p] = rearr(x0)
% Rearranges x0 into coordinates (x) and parameters (p)
global cds nsmds
nap = length(nsmds.ActiveParams);
p = nsmds.P0;
p(nsmds.ActiveParams) = x0((nsmds.nphase+1):end);
x = x0(1:nsmds.nphase);

%HM 060306: Checked the branch switching 

