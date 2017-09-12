clear; clc;
global sys

sys.gui.pausespecial=0;  %Pause at special points 
sys.gui.pausenever=1;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

syshandle=@m2;  %Specify system file

SubFunHandles=feval(syshandle);  %Get function handles from system file
RHShandle=SubFunHandles{2};      %Get function handle for ODE

%Set ODE parameter
k0=.001;k1=10;
k2=0.6;k3=0.2;
k4=1;k5=1;

xinit=[0,0,0,0,0]; %Set ODE initial condition

%Specify ODE function with ODE parameters set
RHS_no_param=@(t,x)RHShandle(t,x,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9); 

%Set ODE integrator parameters.
options=odeset;
options=odeset(options,'RelTol',1e-5);
options=odeset(options,'maxstep',1e-1);

%Integrate until a steady state is found.
[tout xout]=ode45(RHS_no_param,[0,1000],xinit,options);
%figure();plot(xout(:,1))
%
%%%%% Continuation from equilibrium %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set initial condition as the endpoint of integration.  Use
%to bootstrap the continuation.
xinit=xout(size(xout,1),:);

pvec=[k0,k1,k2,k3,k4,k5,k6,k7,k8,k9]';      % Initialize parameter vector

ap=1;

[x0,v0]=init_EP_EP(syshandle, xinit', pvec, ap); %Initialize equilibrium


opt=contset;
opt=contset(opt,'MaxNumPoints',700); %Set numeber of continuation steps
opt=contset(opt,'MaxStepsize',.01);  %Set max step size
opt=contset(opt,'Singularities',1);  %Monitor singularities
opt=contset(opt,'Eigenvalues',1);    %Output eigenvalues 
opt = contset(opt,'InitStepsize',0.01); %Set Initial stepsize

paramid=6;%parameter to vary
ssvalueid=1;%equilibrium level of Al

[x1,v1,s1,h1,f1]=cont(@equilibrium,x0,v0,opt);

% %%
% %%%%% Continuation from equilibrium backward %%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x0,v0]=init_EP_EP(syshandle, xinit', pvec, ap); %Initialize equilibrium
% opt=contset(opt,'Backward',1);
% [x2,v2,s2,h2,f2]=cont(@equilibrium,x0,v0,opt);
% % cpl(x2,v2,s2,[3 1]);
figure();hold on;
plot(x1(paramid,:),x1(ssvalueid,:),'-b');
plot(x1(paramid,cat(1,s1(2:end-1).index)),x1(ssvalueid,cat(1,s1(2:end-1).index)),'ko');

%%
%%%%% Branch swiching and continuation %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chosen=s1(5);
BP=x1(:,chosen.index);
pvec(ap)=BP(paramid);
[x0,vO]=init_BP_EP(syshandle, BP(1:paramid-1), pvec, chosen, 0.01);  
opt=contset(opt,'MaxNumPoints',300);
[x3,v3,s3,h3,f3]=cont(@equilibrium,x0,v0,opt); %Switch branches and continue.

xNearBP=x3(1:paramid-1,10);
pvec(ap)=x3(paramid,10);

[x0,v0]=init_EP_EP(syshandle, xNearBP, pvec, ap); %Initialize equilibrium

opt=contset(opt,'MaxNumPoints',500);
opt=contset(opt,'Backward',1);
[x4,v4,s4,h4,f4]=cont(@equilibrium,x0,v0,opt); %Switch branches and continue.

figure();hold on;
plot(x1(paramid,:),x1(ssvalueid,:),'-b');
plot(x1(paramid,cat(1,s1(2:end-1).index)),x1(ssvalueid,cat(1,s1(2:end-1).index)),'bs');

plot(x3(paramid,:),x3(ssvalueid,:),'-r');
plot(x3(paramid,cat(1,s3(2:end-1).index)),x3(ssvalueid,cat(1,s3(2:end-1).index)),'rs');

plot(x4(paramid,:),x4(ssvalueid,:),'-m');
plot(x4(paramid,cat(1,s4(2:end-1).index)),x4(ssvalueid,cat(1,s4(2:end-1).index)),'ms');

%%
figure();hold on;
plot3(x1(indx,:),x1(indy,:),real(f1(4,:)),'b-')
plot3(x3(indx,:),x3(indy,:),real(f3(4,:)),'r-')
plot3(x4(indx,:),x4(indy,:),real(f4(4,:)),'m-')
%%


x1line = plotStability( x1,f1,paramid-1 );
x3line = plotStability( x3,f3,paramid-1 );
x4line = plotStability( x4,f4,paramid-1 );

figure();hold on;
xline = cat(1,x1line,x3line,x4line);
for i=1:size(xline,1)
    seg=xline{i,1};
    if xline{i,2}
        plot(seg(paramid,:),seg(ssvalueid,:),'r');
    else
        plot(seg(paramid,:),seg(ssvalueid,:),'b');
    end
end
plot(x1(paramid,cat(1,s1(2:end-1).index)),x1(ssvalueid,cat(1,s1(2:end-1).index)),'ko');
text(x1(paramid,cat(1,s1(2:end-1).index))+0.005,x1(ssvalueid,cat(1,s1(2:end-1).index)),cat(1,s1(2:end-1).label));
plot(x3(paramid,cat(1,s3(2:end-1).index)),x3(ssvalueid,cat(1,s3(2:end-1).index)),'ko');
text(x3(paramid,cat(1,s3(2:end-1).index))+0.005,x3(ssvalueid,cat(1,s3(2:end-1).index)),cat(1,s3(2:end-1).label));
plot(x4(paramid,cat(1,s4(2:end-1).index)),x4(ssvalueid,cat(1,s4(2:end-1).index)),'ko');
text(x4(paramid,cat(1,s4(2:end-1).index))+0.005,x4(ssvalueid,cat(1,s4(2:end-1).index)),cat(1,s4(2:end-1).label));
plot(k0,xout(1,1),'c^-');plot(k0,xout(end,1),'cv-');
%%
%Set ODE parameter
k0=.04;
fun=@(x,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9) (k0+k1*x^3/(x^3+k2^3))*(k8/k9-x)-(k3+k4*k5/k6*x)*x;
fun2=@(x) fun(x,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9);

hss(1)=fzero(fun2,1);hss(2)=hss(1);hss(3)=k8/k9; hss(4)=k5/k6*hss(1);hss(5)=hss(4);
delta=0.00001; hss(1)=hss(1)+delta;%hss(2)=hss(2)-delta;
%Specify ODE function with ODE parameters set
RHS_no_param=@(t,x)RHShandle(t,x,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9); 

%Set ODE integrator parameters.
options=odeset;
options=odeset(options,'RelTol',1e-5);
options=odeset(options,'maxstep',1e-1);

%Integrate until a steady state is found.
[tout xout]=ode45(RHS_no_param,[0,500],hss,options);
figure(); hold on; plot(tout,xout(:,1),'r-');plot(tout,xout(:,2),'r:');plot(tout,xout(:,3),'y-');plot(tout,xout(:,4),'b-');plot(tout,xout(:,5),'b:');
%%
i=3;
k0=params(i,1);k1=params(i,2);k2=params(i,3);k3=params(i,4);k4=params(i,5);k5=params(i,6);

A=linspace(0,10,1000);
F=((k0+k1*A.^3/(A.^3+1))*k4/k5./A-k2)/k3;

zx = zci(F-A);% Approximate Zero-Crossing Indices
figure();hold on;plot(A,F);plot(A,A);plot(A(zx),F(zx),'kx');axis([0 inf -10 10]);
%%
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);% Returns Zero-Crossing Indices Of Argument Vector
K0=10;K1=10;
K2=10;K3=10;
K4=1;K5=1;
A=linspace(0,10,1000);
params=[];N=0;
while N<10
    k0=K0*rand();k1=K1*rand();k2=K2*rand();k3=K2*rand();k4=K4;k5=K5;
    F=((k0+k1*A.^3/(A.^3+1))*k4/k5./A-k2)/k3;%size(zci(F-A),1)
    if size(zci(F-A),1)>2
        N=N+1; 
        params(N,:)=[k0,k1,k2,k3,k4,k5];
        [k0,k1,k2,k3,k4,k5]
    end
end