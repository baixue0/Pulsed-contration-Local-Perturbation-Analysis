%% Test R1 Normal Form
%     clc
%     clear all
%     close all
%     init
%     %% Fixed point continuation
%     opt = contset;
%     p=[-.1,.1,1,.2];ap=1;
%     [x0,v0]=init_FPm_FPm(@NFR1,[-.4;0], p, ap,1);
%     opt=contset(opt,'Backward',0);
%     opt=contset(opt,'Singularities',1);
%     [x1,v1,s1,h1,f1]=cont(@fixedpointmap,x0,[],opt);
%     p=[-.1,-.1,1,.2];ap=2;
%     [x0,v0]=init_FPm_FPm(@NFR1,[-0.3;0.0], p, ap,1);
%     opt=contset(opt,'Backward',1);
%     [x2,v2,s2,h2,f2]=cont(@fixedpointmap,x0,[],opt);
%     
%     %% LP and NS continuation
%     x0=x1(1:2,s1(3).index);p0=p;p0(1)=x1(3,s1(3).index);ap=[1 2]; 
%     [x0,v0]=init_LPm_LPm(@NFR1,x0, p0, ap,1);
%     opt=contset(opt,'Backward',0);
%     [xlp1,vlp1,slp1,hlp1,flp1]=cont(@limitpointmap,x0,[],opt);
%     opt=contset(opt,'Backward',1);
%     [xlp2,vlp2,slp2,hlp2,flp2]=cont(@limitpointmap,x0,[],opt);
% 
%     x0=x2(1:2,s2(2).index);p0=p;p0(2)=x2(3,s2(2).index);ap=[1 2]; 
%     [x0,v0]=init_NSm_NSm(@NFR1,x0, p0', ap,1);
%     opt=contset(opt,'IgnoreSingularity',[2 3 5]);
%     opt=contset(opt,'Backward',0);
%     [xns1,vns1,sns1,hns1,fns1]=cont(@neimarksackermap,x0,[],opt);
%     opt=contset(opt,'Backward',1);
%     [xns2,vns2,sns2,hns2,fns2]=cont(@neimarksackermap,x0,[],opt);
%     
%     %%Plot these results;
%     figure(1);
%     cpl(xlp1,vlp1,slp1,[3  4]); cpl(xlp2,vlp2,slp2,[3  4]);
%     cpl(xns1,vns1,sns1,[3  4]); cpl(xns2,vns2,sns2,[3  4]);
% 

    clc
    global opt
    opt=contset;
    p0=[-0.1,-.383031,1,.2];
%     p0=[-0.2,-.53125,1,.2];   
%     p0=[-0.5,-.815,1,.2];
    x0=[sqrt(-p0(1)/p0(3));0];

    optM=init_FPm_1DMan(@NFR1,x0,p0,1);
    %% Uman
    optM.function='Umanifold';  
    % I varied parameters in different continuations and in order to check 
    % the results
    optM.eps=1e-5;
    optM.nmax=40000;
    optM.deltaMax=0.0001;
    optM.deltak=1e-5;
    optM.deltaAlphaMax=0.00001;
    optM.deltaAlphaMin=0.000001;
    optM.Niterations=1;
    
    optM.distanceInit=-optM.distanceInit;
    [a,l]=contman(optM);toc
    optM.distanceInit=-optM.distanceInit;
    [a1,l1]=contman(optM);toc
    %% Sman
    optM=init_FPm_1DMan(@NFR1,x0,p0,1);
    optM.eps=1e-8;
    optM.Niterations=1;
    optM.deltaMax=0.0001;
    optM.deltak=1e-5;
    optM.deltaAlphaMin=0.00001;
    optM.function='Smanifold';    
    % n.b. some parameters are used only here, ex: NtwMax
    optM.distanceInit=-optM.distanceInit;
    optM.nmax=1000;
    [b,l]=contman(optM);toc
    line(b(1,:),b(2,:),'color','b')
    optM.distanceInit=-optM.distanceInit;
    [b1,l]=contman(optM);toc
    %% plot
    figure(2);clf(2);
    ylim([-2,2]);
    xlim([-2,2]);
    line(x0(1),x0(2),'marker','.','color','k')
    line(a(1,:),a(2,:),'color','r')
    line(a1(1,:),a1(2,:),'color','r')
    line(b(1,:),b(2,:),'color','b')  
    line(b1(1,:),b1(2,:),'color','b')
%     %% homocline curves
%     hom=findintersections(a,b)
%     %% plot some homoclinic curves
%     % the primaries (the 2 longest curves)
%     c=hom{4};
%     d=hom{end};
%     line(c(1,:),c(2,:),'color','g','marker','.')
%     line(d(1,:),d(2,:),'color','m','marker','.')
%     % some secondaries
%     e=hom{12};
%     f=hom{13};
%     g=hom{14};
%     h=hom{15};
%     line(e(1,:),e(2,:),'color','g','marker','.')
%     line(f(1,:),f(2,:),'color','m','marker','.')
%     line(g(1,:),g(2,:),'color','c','marker','.')
%     line(h(1,:),h(2,:),'color','y','marker','.')
%     %% continuation of the primary curve
%     opt=contset;
%     opt=contset(opt,'MaxNumPoints',100);
%     opt=contset(opt,'Singularities',1);
%     opt=contset(opt,'MinStepsize',1e-8);
%     [x0hom,v0hom]=init_Hom_Hom(@Ghmap,[x0,c],p0,2,1);
%     [xhomc1,vhomc1,shomc1,hhomc1,fhomc1]=cont(@homoclinic,x0hom,[],opt);
%     opt=contset(opt,'Backward',1);
%     [xhomc2,vhomc2,shomc2,hhomc2,fhomc2]=cont(@homoclinic,x0hom,[],opt);
%     
%     figure
%     cpl(xhomc1,vhomc1,shomc1,[size(xhomc1,1),14])
%     cpl(xhomc2,vhomc2,shomc2,[size(xhomc2,1),14])
%     
% func_handles = feval(@Ghmap);
% [V,D]=eig(cjac(@Ghmap,func_handles{3},x0,n2c(p0)));
% lu=D(2,2);ls=D(1,1);
% Nu_ext=15;Ns_ext=15;
% cextend=[(c(:,end)-x0)*lu.^(-Nu_ext:-1)+repmat(x0,1,Nu_ext) ...
%           c (c(:,end)-x0)*ls.^(1:Ns_ext)+repmat(x0,1,Ns_ext)];
% opt=contset(opt,'Singularities',1);
% opt = contset(opt,'Backward',0);
% opt=contset(opt,'MaxNumPoints',40);
% [x0hom,v0hom]=init_Hom_Hom(@Ghmap,[x0,cextend],p0,2,1);
% [xhomc1,vhomc1,shomc1,hhomc1,fhomc1]=cont(@homoclinic,x0hom,[],opt);
% 
% opt=contset(opt,'Singularities',0);
% opt = contset(opt,'Backward',1);
% opt=contset(opt,'MaxNumPoints',10);
% p1=p0;p1(2)=xhomc1(end,shomc1(2).index);
% [x0lphom,v0lphom]=init_HomT_HomT(@Ghmap,xhomc1(1:end-1,shomc1(2).index),2,1,1,p1,[1 2],1);
% [xhomT,vhomT,shomT,homT,fhomT]=cont(@homoclinicT,x0lphom,[],opt);
% cpl(xhomT,vhomT,shomT,[size(xhomT,1),size(xhomT,1)-1]); %plot in (beta,alpha) plane
