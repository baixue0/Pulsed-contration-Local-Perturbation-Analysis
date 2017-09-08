opt = contset;
opt = contset(opt,'Singularities',1);
opt = contset(opt,'Multipliers',1);
%opt = contset(opt,'MaxNumPoints',1250);
p=[0;-2.5;-.1;0];ap=1;
[x0,v0]=init_FPm_FPm(@Ghmap,[-2;-2], p, ap,1);
[x1,v1,s1,h1,f1]=cont(@fixedpointmap,x0,[],opt,1);

%Limit point; first iterate
xlp=x1(1:2,s1(3).index);plp=p;plp(ap) = x1(end,s1(3).index);alp = [1 2];
opt = contset(opt,'Backward',1);
[x0,v0]=init_LPm_LPm(@Ghmap,xlp,plp,alp,1);
[x2,v2,s2,h2,f2]=cont(@limitpointmap,x0,[],opt,1);
cpl(x2,v2,s2,[4 3]);

%Period doubling; first iterate
xpd=x1(1:2,s1(2).index);ppd=p;ppd(ap) = x1(end,s1(2).index);apd = [1 2];
opt = contset(opt,'Backward',0);
[x0,v0]=init_PDm_PDm(@Ghmap,xpd,ppd,apd,1);
[x3,v3,s3,h3,f3]=cont(@perioddoublingmap,x0,[],opt,1);
cpl(x3,v3,s3,[4 3]);

%Neimark-Sacker; first iterate
xns=x2(1:2,s2(3).index);pns=p;pns(alp) = x2(3:4,s2(3).index);aans = [1 2];
opt = contset(opt,'Backward',0);
opt=contset(opt,'IgnoreSingularity',[5]);
[x0,v0]=init_NSm_NSm(@Ghmap,xns,pns,aans,1);
[x4,v4,s4,h4,f4]=cont(@neimarksackermap,x0,[],opt,1);
[x4,v4,s4,h4,f4]=cont(@neimarksackermap,x0,[],opt,1);
cpl(x4,v4,s4,[4 3]);

%Neimark-Sacker; second iterate; starting from LPPD
xns2=x2(1:2,s2(2).index);pns2=plp;pns2(alp) = x2(3:4,s2(2).index);ans2 = [1 2];
opt = contset(opt,'Backward',1);
opt = contset(opt,'MaxStepSize', .05);
opt = contset(opt,'MaxNumPoints',150);
opt=contset(opt,'IgnoreSingularity',[5]);
[x0,v0]=init_LPPD_NS2m(@Ghmap,1e-3,xns2,pns2,ans2,1);
[x5,v5,s5,h5,f5]=cont(@neimarksackermap,x0,[],opt,2);
cpl(x5,v5,s5,[4 3]);

%Neimark-Sacker; second iterate; starting from R2
xns2=x4(1:2,s4(5).index);pns2=pns;pns2(aans) = x4(3:4,s4(5).index);ans2 = [1 2];
opt = contset(opt,'Backward',1);
opt = contset(opt,'MaxNumPoints',150);
opt=contset(opt,'IgnoreSingularity',[5]);
[x0,v0]=init_R2_NS2m(@Ghmap,1e-3,xns2,pns2,ans2,1);
[x6,v6,s6,h6,f6]=cont(@neimarksackermap,x0,[],opt,2);
cpl(x6,v6,s6,[4 3]);

%Period doubling; second iterate; starting from R2
xpd2=x6(1:2,s6(4).index);ppd2=pns2;ppd2(ans2) = x6(3:4,s6(4).index);apd2 = [1 2];
opt = contset(opt,'Backward',1);
opt = contset(opt,'MaxNumPoints',150);
[x0,v0]=init_PDm_PDm(@Ghmap,xpd2,ppd2,apd2,2);
[x7,v7,s7,h7,f7]=cont(@perioddoublingmap,x0,[],opt,2);
cpl(x7,v7,s7,[4 3]);
opt = contset(opt,'Backward',0);
[x8,v8,s8,h8,f8]=cont(@perioddoublingmap,x0,[],opt,2);
cpl(x8,v8,s8,[4 3]);

%Limit Point and Neimark-Sacker; fourth iterate; starting from R4
xlp4=x4(1:2,s4(3).index);plp4=pns;plp4(aans) = x4(3:4,s4(3).index);alp4 = [1 2];
[x0,v0]=init_R4_LP4m1(@Ghmap,1e-3,xlp4,plp4,alp4,1);
opt = contset(opt,'MaxNumPoints',250);
opt = contset(opt,'Backward',0);
[x9,v9,s9,h9,f9]=cont(@limitpointmap,x0,[],opt,4);
cpl(x9,v9,s9,[4 3]);
[x0,v0]=init_R4_LP4m2(@Ghmap,1e-3,xlp4,plp4,alp4,1);
opt = contset(opt,'Backward',0);
[x10,v10,s10,h10,f10]=cont(@limitpointmap,x0,[],opt,4);
cpl(x10,v10,s10,[4 3]);
[x0,v0]=init_R4_NS4m(@Ghmap,1e-3,xlp4,plp4,alp4,1);%%%%WRONG
opt = contset(opt,'Backward',1);
[x11,v11,s11,h11,f11]=cont(@neimarksackermap,x0,[],opt,4);
cpl(x11,v11,s11,[4 3]);
%Period doubling; fourth iterate; starting from R2
xpd4=x11(1:2,s11(4).index);ppd4=pns;ppd4(aans) = x11(3:4,s11(4).index);apd4 =[1 2];
[x0,v0]=init_PDm_PDm(@Ghmap,xpd4,ppd4,apd4,4);
opt = contset(opt,'Backward',0);
[x12,v12,s12,h12,f12]=cont(@perioddoublingmap,x0,[],opt,4);
opt = contset(opt,'Backward',1);
[x13,v13,s13,h13,f13]=cont(@perioddoublingmap,x0,[],opt,4);
cpl(x12,v12,s12,[4 3]);cpl(x13,v13,s13,[4 3]);

%Neutral Saddle; third iterate; starting from R3
xns3=x4(1:2,s4(4).index);pns3=pns;pns3(aans) = x4(3:4,s4(4).index);ans3 = [1 2];
[x0,v0]=init_R3_NS3m(@Ghmap,1e-4,xns3,pns3,ans3,1);
[x14,v14,s14,h14,f14]=cont(@neimarksackermap,x0,[],opt,3);
cpl(x14,v14,s14,[4 3]);

%Limit Point; third iterate; starting from R1
xlp3=x14(1:2,s14(2).index);plp3=pns;plp3(aans) = x14(3:4,s14(2).index);alp3 =[1 2];
[x0,v0]=init_LPm_LPm(@Ghmap,xlp3,plp3,alp3,3);
[x15,v15,s15,h15,f15]=cont(@limitpointmap,x0,[],opt,3);
opt = contset(opt,'Backward',0);
[x16,v16,s16,h16,f16]=cont(@limitpointmap,x0,[],opt,3);
cpl(x15,v15,s15,[4 3]);cpl(x16,v16,s16,[4 3]);
%Period doubling; third iterate; starting from R2
xpd3=x14(1:2,s14(5).index);ppd3=pns;ppd3(aans) = x14(3:4,s14(5).index);apd3 =[1 2];
[x0,v0]=init_PDm_PDm(@Ghmap,xpd3,ppd3,apd3,3);
[x17,v17,s17,h17,f17]=cont(@perioddoublingmap,x0,[],opt,3);
opt = contset(opt,'Backward',1);
[x18,v18,s18,h18,f18]=cont(@perioddoublingmap,x0,[],opt,3);
cpl(x17,v17,s17,[4 3]);cpl(x18,v18,s18,[4 3]);


%Degenerate flip: Switching twice
opt = contset(opt,'Backward',0);
xlp8=x13(1:2,s13(2).index);plp8=pns;plp8(aans) = x13(3:4,s13(2).index);alp8 =[1 2];
[x0,v0]=init_GPD_LP2m(@Ghmap,1e-4,xlp8,plp8,alp8,4);
[x20,v20,s20,h20,f20]=cont(@limitpointmap,x0,[],opt,8);
opt = contset(opt,'MaxNumPoints',800);
xlp8=x13(1:2,s13(3).index);plp8=pns;plp8(aans) = x13(3:4,s13(3).index);alp8 =[1 2];
[x0,v0]=init_GPD_LP2m(@Ghmap,1e-4,xlp8,plp8,alp8,4);
[x21,v21,s21,h21,f21]=cont(@limitpointmap,x0,[],opt,8);
cpl(x20,v20,s20,[4 3]);cpl(x21,v21,s21,[4 3]);


xns8=x20(1:2,s20(3).index);pns8=pns;pns8(aans) = x20(3:4,s20(3).index);ans8 =[1 2];
[x0,v0]=init_NSm_NSm(@Ghmap,xns8,pns8,ans8,8);
%[x22,v22,s22,h22,f22]=cont(@neimarksackermap,x0,[],opt,8);
%cpl(x22,v22,s22,[4 3]);

%The large figure
figure(2);axis([-1.5 1.5 -1 5]);
cpl(x2,v2,s2,[4,3]);cpl(x3,v3,s3,[4,3]);
cpl(x4(:,1:s4(5).index),v4(:,1:s4(5).index),s4(1:5),[4 3]);
cpl(x5(:,1:s5(4).index),v5(:,1:s5(4).index),s5(1:4),[4 3]);
cpl(x6(:,1:s6(4).index),v6(:,1:s6(4).index),s6(1:4),[4 3]);
cpl(x7,v7,s7,[4 3]);cpl(x8,v8,s8,[4 3]);
cpl(x9,v9,s9,[4 3]);cpl(x10,v10,s10,[4 3]);
cpl(x11(:,1:s11(4).index),v11(:,1:s11(4).index),s11(1:4),[4 3]);
cpl(x12,v12,s12,[4 3]);cpl(x13,v13,s13,[4 3]);
cpl(x14(:,1:s14(5).index),v11(:,1:s14(5).index),s14(1:5),[4 3]);
cpl(x15,v15,s15,[4 3]);cpl(x16,v16,s16,[4 3]);cpl(x17,v17,s17,[4 3]);cpl(x18,v18,s18,[4 3]);
cpl(x20,v20,s20,[4 3]);cpl(x21,v21,s21,[4 3]);

%Zoom near R4
figure(3);axis([.98 1.02 -.04 .06]);
cpl(x4(:,1:s4(5).index),v4(:,1:s4(5).index),s4(1:5),[4 3]);
cpl(x9,v9,s9,[4 3]);cpl(x10,v10,s10,[4 3]);
cpl(x11(:,1:s11(4).index),v11(:,1:s11(4).index),s11(1:4),[4 3]);
cpl(x12,v12,s12,[4 3]);cpl(x13,v13,s13,[4 3]);

steps=[1e-7;1e-6;1e-5;5e-5;1e-4;5e-4;1e-3;2e-3;5e-3;1e-2;5e-2;.1;.5];
predict_R4_1=zeros(4,12);predict_R4_1(:,1)=[xlp4;plp4(aans)];
predict_R4_2=zeros(4,12);predict_R4_2(:,1)=[xlp4;plp4(aans)];
predict_R4_3=zeros(4,12);predict_R4_3(:,1)=[xlp4;plp4(aans)];
global cds
for i=1:11
    cds.options.AutDerivative=1;cds.options.AutDerivativeIte=24;
    [x0,v0]=init_R4_LP4m1(@Ghmap,steps(i),xlp4,plp4,alp4,1);
    predict_R4_1(:,i+1)=x0;
    cds.options.AutDerivative=1;cds.options.AutDerivativeIte=24;
    [x0,v0]=init_R4_LP4m2(@Ghmap,steps(i),xlp4,plp4,alp4,1);
    predict_R4_2(:,i+1)=x0;
    cds.options.AutDerivative=1;cds.options.AutDerivativeIte=24;
    [x0,v0]=init_R4_NS4m(@Ghmap,steps(i),xlp4,plp4,alp4,1);
    predict_R4_3(:,i+1)=x0(1:4);
end
hold on;plot(predict_R4_1(4,:),predict_R4_1(3,:),'linestyle','-.')
hold on;plot(predict_R4_2(4,:),predict_R4_2(3,:),'linestyle','-.')
hold on;plot(predict_R4_3(4,:),predict_R4_3(3,:),'linestyle','-.')
figure(13);%check phase too
cpl(x9,v9,s9,[1 2]);cpl(x10,v10,s10,[1 2]);cpl(x11,v11,s11,[1 2]);
hold on;plot(predict_R4_1(1,:),predict_R4_1(2,:),'linestyle','-.')
hold on;plot(predict_R4_2(1,:),predict_R4_2(2,:),'linestyle','-.')
hold on;plot(predict_R4_3(1,:),predict_R4_3(2,:),'linestyle','-.')

%Zoom near R3
figure(4);axis([.9 1.05 .9 1.3]);
cpl(x4(:,1:s4(5).index),v4(:,1:s4(5).index),s4(1:5),[4 3]);
cpl(x14(:,1:s14(5).index),v11(:,1:s14(5).index),s14(1:5),[4 3]);
cpl(x15,v15,s15,[4 3]);cpl(x16,v16,s16,[4 3]);cpl(x17,v17,s17,[4 3]);cpl(x18,v18,s18,[4 3]);

steps=[1e-5;5e-5;1e-4;5e-4;1e-3;2e-3;5e-3;8e-3;.01;.03;.05;.1;.15;.25;.4];
predict_R3_1=zeros(4,16);predict_R3_1(:,1)=[xns3;pns3(aans)];
predict_R3_2=zeros(4,16);predict_R3_2(:,1)=[xns3;pns3(aans)];
for i=1:15
    cds.options.AutDerivative=1;cds.options.AutDerivativeIte=24;
    [x0,v0]=init_R3_NS3m(@Ghmap,-steps(i),xns3,pns3,ans3,1);
    predict_R3_1(:,i+1)=x0(1:4);
    cds.options.AutDerivative=1;cds.options.AutDerivativeIte=24;
    [x0,v0]=init_R3_NS3m(@Ghmap,steps(i),xns3,pns3,ans3,1);
    predict_R3_2(:,i+1)=x0(1:4);
end
hold on;plot(predict_R3_1(4,:),predict_R3_1(3,:),'linestyle','-.')
hold on;plot(predict_R3_2(4,:),predict_R3_2(3,:),'linestyle','-.')

%Zoom near DPD
figure(5);axis([-1 -.8 2.6 3]);
cpl(x9,v9,s9,[4 3]);cpl(x10,v10,s10,[4 3]);
cpl(x12,v12,s12,[4 3]);cpl(x13,v13,s13,[4 3]);
cpl(x20,v20,s20,[4 3]);cpl(x21,v21,s21,[4 3]);

steps=[1e-5;5e-5;1e-4;5e-4;1e-3;5e-3;.01;.04;.07;.1;.15;.2;.3;.4];
xlp81=x13(1:2,s13(2).index);plp81=pns;plp81(ans3) = x13(3:4,s13(2).index);
xlp82=x13(1:2,s13(3).index);plp82=pns;plp82(ans3) = x13(3:4,s13(3).index);
predict_LP_1=zeros(4,15);predict_LP_1(:,1)=[xlp81;plp81(aans)];
predict_LP_2=zeros(4,15);predict_LP_2(:,1)=[xlp82;plp82(aans)];
for i=1:14
    cds.options.AutDerivative=1;cds.options.AutDerivativeIte=24;
    [x0,v0]=init_GPD_LP2m(@Ghmap,steps(i),xlp81,plp81,ans3,4);
    predict_LP_1(:,i+1)=x0(1:4);
    cds.options.AutDerivative=1;cds.options.AutDerivativeIte=24;
    [x0,v0]=init_GPD_LP2m(@Ghmap,steps(i),xlp82,plp82,ans3,4);
    predict_LP_2(:,i+1)=x0(1:4);
end
hold on;plot(predict_LP_1(4,:),predict_LP_1(3,:),'linestyle','-.')
hold on;plot(predict_LP_2(4,:),predict_LP_2(3,:),'linestyle','-.')
