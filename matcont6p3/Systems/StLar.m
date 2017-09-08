function out = StLar
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,km7,k8)
dydt=[-k1*kmrgd(1)*kmrgd(2)*kmrgd(3)-k3*kmrgd(1)*kmrgd(2)*kmrgd(4)+k7-km7*kmrgd(1);;
-k1*kmrgd(1)*kmrgd(2)*kmrgd(3)-k3*kmrgd(1)*kmrgd(2)*kmrgd(4)+k8;;
k1*kmrgd(1)*kmrgd(2)*kmrgd(3)-2*k2*kmrgd(3)^2+2*k3*kmrgd(1)*kmrgd(2)*kmrgd(4)-k4*kmrgd(3)+k6;;
-k3*kmrgd(1)*kmrgd(2)*kmrgd(4)+2*k2*kmrgd(3)^2-k5*kmrgd(4);;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(StLar);
y0=[0,0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,km7,k8)
jac=[ - km7 - kmrgd(2)*kmrgd(3)*k1 - kmrgd(2)*kmrgd(4)*k3 , - kmrgd(1)*kmrgd(3)*k1 - kmrgd(1)*kmrgd(4)*k3 , -kmrgd(1)*kmrgd(2)*k1 , -kmrgd(1)*kmrgd(2)*k3 ; - kmrgd(2)*kmrgd(3)*k1 - kmrgd(2)*kmrgd(4)*k3 , - kmrgd(1)*kmrgd(3)*k1 - kmrgd(1)*kmrgd(4)*k3 , -kmrgd(1)*kmrgd(2)*k1 , -kmrgd(1)*kmrgd(2)*k3 ; kmrgd(2)*kmrgd(3)*k1 + 2*kmrgd(2)*kmrgd(4)*k3 , kmrgd(1)*kmrgd(3)*k1 + 2*kmrgd(1)*kmrgd(4)*k3 , kmrgd(1)*kmrgd(2)*k1 - 4*kmrgd(3)*k2 - k4 , 2*kmrgd(1)*kmrgd(2)*k3 ; -kmrgd(2)*kmrgd(4)*k3 , -kmrgd(1)*kmrgd(4)*k3 , 4*kmrgd(3)*k2 , - k5 - kmrgd(1)*kmrgd(2)*k3 ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,km7,k8)
jacp=[ -kmrgd(1)*kmrgd(2)*kmrgd(3) , 0 , -kmrgd(1)*kmrgd(2)*kmrgd(4) , 0 , 0 , 0 , 1 , -kmrgd(1) , 0 ; -kmrgd(1)*kmrgd(2)*kmrgd(3) , 0 , -kmrgd(1)*kmrgd(2)*kmrgd(4) , 0 , 0 , 0 , 0 , 0 , 1 ; kmrgd(1)*kmrgd(2)*kmrgd(3) , -2*kmrgd(3)^2 , 2*kmrgd(1)*kmrgd(2)*kmrgd(4) , -kmrgd(3) , 0 , 1 , 0 , 0 , 0 ; 0 , 2*kmrgd(3)^2 , -kmrgd(1)*kmrgd(2)*kmrgd(4) , 0 , -kmrgd(4) , 0 , 0 , 0 , 0 ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,km7,k8)
hess1=[ 0 , - kmrgd(3)*k1 - kmrgd(4)*k3 , -kmrgd(2)*k1 , -kmrgd(2)*k3 ; 0 , - kmrgd(3)*k1 - kmrgd(4)*k3 , -kmrgd(2)*k1 , -kmrgd(2)*k3 ; 0 , kmrgd(3)*k1 + 2*kmrgd(4)*k3 , kmrgd(2)*k1 , 2*kmrgd(2)*k3 ; 0 , -kmrgd(4)*k3 , 0 , -kmrgd(2)*k3 ];
hess2=[ - kmrgd(3)*k1 - kmrgd(4)*k3 , 0 , -kmrgd(1)*k1 , -kmrgd(1)*k3 ; - kmrgd(3)*k1 - kmrgd(4)*k3 , 0 , -kmrgd(1)*k1 , -kmrgd(1)*k3 ; kmrgd(3)*k1 + 2*kmrgd(4)*k3 , 0 , kmrgd(1)*k1 , 2*kmrgd(1)*k3 ; -kmrgd(4)*k3 , 0 , 0 , -kmrgd(1)*k3 ];
hess3=[ -kmrgd(2)*k1 , -kmrgd(1)*k1 , 0 , 0 ; -kmrgd(2)*k1 , -kmrgd(1)*k1 , 0 , 0 ; kmrgd(2)*k1 , kmrgd(1)*k1 , -4*k2 , 0 ; 0 , 0 , 4*k2 , 0 ];
hess4=[ -kmrgd(2)*k3 , -kmrgd(1)*k3 , 0 , 0 ; -kmrgd(2)*k3 , -kmrgd(1)*k3 , 0 , 0 ; 2*kmrgd(2)*k3 , 2*kmrgd(1)*k3 , 0 , 0 ; -kmrgd(2)*k3 , -kmrgd(1)*k3 , 0 , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
hess(:,:,4) =hess4;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,km7,k8)
hessp1=[ -kmrgd(2)*kmrgd(3) , -kmrgd(1)*kmrgd(3) , -kmrgd(1)*kmrgd(2) , 0 ; -kmrgd(2)*kmrgd(3) , -kmrgd(1)*kmrgd(3) , -kmrgd(1)*kmrgd(2) , 0 ; kmrgd(2)*kmrgd(3) , kmrgd(1)*kmrgd(3) , kmrgd(1)*kmrgd(2) , 0 ; 0 , 0 , 0 , 0 ];
hessp2=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , -4*kmrgd(3) , 0 ; 0 , 0 , 4*kmrgd(3) , 0 ];
hessp3=[ -kmrgd(2)*kmrgd(4) , -kmrgd(1)*kmrgd(4) , 0 , -kmrgd(1)*kmrgd(2) ; -kmrgd(2)*kmrgd(4) , -kmrgd(1)*kmrgd(4) , 0 , -kmrgd(1)*kmrgd(2) ; 2*kmrgd(2)*kmrgd(4) , 2*kmrgd(1)*kmrgd(4) , 0 , 2*kmrgd(1)*kmrgd(2) ; -kmrgd(2)*kmrgd(4) , -kmrgd(1)*kmrgd(4) , 0 , -kmrgd(1)*kmrgd(2) ];
hessp4=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , -1 , 0 ; 0 , 0 , 0 , 0 ];
hessp5=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , -1 ];
hessp6=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp7=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp8=[ -1 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp9=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
hessp(:,:,6) =hessp6;
hessp(:,:,7) =hessp7;
hessp(:,:,8) =hessp8;
hessp(:,:,9) =hessp9;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,km7,k8)
tens31=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens32=[ 0 , 0 , -k1 , -k3 ; 0 , 0 , -k1 , -k3 ; 0 , 0 , k1 , 2*k3 ; 0 , 0 , 0 , -k3 ];
tens33=[ 0 , -k1 , 0 , 0 ; 0 , -k1 , 0 , 0 ; 0 , k1 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens34=[ 0 , -k3 , 0 , 0 ; 0 , -k3 , 0 , 0 ; 0 , 2*k3 , 0 , 0 ; 0 , -k3 , 0 , 0 ];
tens35=[ 0 , 0 , -k1 , -k3 ; 0 , 0 , -k1 , -k3 ; 0 , 0 , k1 , 2*k3 ; 0 , 0 , 0 , -k3 ];
tens36=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens37=[ -k1 , 0 , 0 , 0 ; -k1 , 0 , 0 , 0 ; k1 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens38=[ -k3 , 0 , 0 , 0 ; -k3 , 0 , 0 , 0 ; 2*k3 , 0 , 0 , 0 ; -k3 , 0 , 0 , 0 ];
tens39=[ 0 , -k1 , 0 , 0 ; 0 , -k1 , 0 , 0 ; 0 , k1 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens310=[ -k1 , 0 , 0 , 0 ; -k1 , 0 , 0 , 0 ; k1 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens311=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens312=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens313=[ 0 , -k3 , 0 , 0 ; 0 , -k3 , 0 , 0 ; 0 , 2*k3 , 0 , 0 ; 0 , -k3 , 0 , 0 ];
tens314=[ -k3 , 0 , 0 , 0 ; -k3 , 0 , 0 , 0 ; 2*k3 , 0 , 0 , 0 ; -k3 , 0 , 0 , 0 ];
tens315=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens316=[ 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,1,3) =tens33;
tens3(:,:,1,4) =tens34;
tens3(:,:,2,1) =tens35;
tens3(:,:,2,2) =tens36;
tens3(:,:,2,3) =tens37;
tens3(:,:,2,4) =tens38;
tens3(:,:,3,1) =tens39;
tens3(:,:,3,2) =tens310;
tens3(:,:,3,3) =tens311;
tens3(:,:,3,4) =tens312;
tens3(:,:,4,1) =tens313;
tens3(:,:,4,2) =tens314;
tens3(:,:,4,3) =tens315;
tens3(:,:,4,4) =tens316;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,km7,k8)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,k1,k2,k3,k4,k5,k6,k7,km7,k8)
