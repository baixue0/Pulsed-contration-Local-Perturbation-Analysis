function out = EcoMod
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
function dydt = fun_eval(t,kmrgd,RR,AA,BB,CC,DD)
dydt=[RR*kmrgd(1)*(1-kmrgd(1)) - kmrgd(1)*kmrgd(2)/(kmrgd(1) + AA);
-CC*kmrgd(2) + kmrgd(1)*kmrgd(2)/(kmrgd(1)+AA) - DD*kmrgd(2)*kmrgd(2)/(kmrgd(2)*kmrgd(2) + BB*BB);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(EcoMod);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,RR,AA,BB,CC,DD)
jac=[ (kmrgd(1)*kmrgd(2))/(AA + kmrgd(1))^2 - kmrgd(2)/(AA + kmrgd(1)) - RR*(kmrgd(1) - 1) - RR*kmrgd(1) , -kmrgd(1)/(AA + kmrgd(1)) ; kmrgd(2)/(AA + kmrgd(1)) - (kmrgd(1)*kmrgd(2))/(AA + kmrgd(1))^2 , kmrgd(1)/(AA + kmrgd(1)) - CC - (2*DD*kmrgd(2))/(BB^2 + kmrgd(2)^2) + (2*DD*kmrgd(2)^3)/(BB^2 + kmrgd(2)^2)^2 ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,RR,AA,BB,CC,DD)
jacp=[ -kmrgd(1)*(kmrgd(1) - 1) , (kmrgd(1)*kmrgd(2))/(AA + kmrgd(1))^2 , 0 , 0 , 0 ; 0 , -(kmrgd(1)*kmrgd(2))/(AA + kmrgd(1))^2 , (2*BB*DD*kmrgd(2)^2)/(BB^2 + kmrgd(2)^2)^2 , -kmrgd(2) , -kmrgd(2)^2/(BB^2 + kmrgd(2)^2) ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,RR,AA,BB,CC,DD)
hess1=[ (2*kmrgd(2))/(AA + kmrgd(1))^2 - 2*RR - (2*kmrgd(1)*kmrgd(2))/(AA + kmrgd(1))^3 , kmrgd(1)/(AA + kmrgd(1))^2 - 1/(AA + kmrgd(1)) ; (2*kmrgd(1)*kmrgd(2))/(AA + kmrgd(1))^3 - (2*kmrgd(2))/(AA + kmrgd(1))^2 , 1/(AA + kmrgd(1)) - kmrgd(1)/(AA + kmrgd(1))^2 ];
hess2=[ kmrgd(1)/(AA + kmrgd(1))^2 - 1/(AA + kmrgd(1)) , 0 ; 1/(AA + kmrgd(1)) - kmrgd(1)/(AA + kmrgd(1))^2 , (10*DD*kmrgd(2)^2)/(BB^2 + kmrgd(2)^2)^2 - (2*DD)/(BB^2 + kmrgd(2)^2) - (8*DD*kmrgd(2)^4)/(BB^2 + kmrgd(2)^2)^3 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,RR,AA,BB,CC,DD)
hessp1=[ 1 - 2*kmrgd(1) , 0 ; 0 , 0 ];
hessp2=[ kmrgd(2)/(AA + kmrgd(1))^2 - (2*kmrgd(1)*kmrgd(2))/(AA + kmrgd(1))^3 , kmrgd(1)/(AA + kmrgd(1))^2 ; (2*kmrgd(1)*kmrgd(2))/(AA + kmrgd(1))^3 - kmrgd(2)/(AA + kmrgd(1))^2 , -kmrgd(1)/(AA + kmrgd(1))^2 ];
hessp3=[ 0 , 0 ; 0 , (4*BB*DD*kmrgd(2))/(BB^2 + kmrgd(2)^2)^2 - (8*BB*DD*kmrgd(2)^3)/(BB^2 + kmrgd(2)^2)^3 ];
hessp4=[ 0 , 0 ; 0 , -1 ];
hessp5=[ 0 , 0 ; 0 , (2*kmrgd(2)^3)/(BB^2 + kmrgd(2)^2)^2 - (2*kmrgd(2))/(BB^2 + kmrgd(2)^2) ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
hessp(:,:,5) =hessp5;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,RR,AA,BB,CC,DD)
tens31=[ (6*kmrgd(1)*kmrgd(2))/(AA + kmrgd(1))^4 - (6*kmrgd(2))/(AA + kmrgd(1))^3 , 2/(AA + kmrgd(1))^2 - (2*kmrgd(1))/(AA + kmrgd(1))^3 ; (6*kmrgd(2))/(AA + kmrgd(1))^3 - (6*kmrgd(1)*kmrgd(2))/(AA + kmrgd(1))^4 , (2*kmrgd(1))/(AA + kmrgd(1))^3 - 2/(AA + kmrgd(1))^2 ];
tens32=[ 2/(AA + kmrgd(1))^2 - (2*kmrgd(1))/(AA + kmrgd(1))^3 , 0 ; (2*kmrgd(1))/(AA + kmrgd(1))^3 - 2/(AA + kmrgd(1))^2 , 0 ];
tens33=[ 2/(AA + kmrgd(1))^2 - (2*kmrgd(1))/(AA + kmrgd(1))^3 , 0 ; (2*kmrgd(1))/(AA + kmrgd(1))^3 - 2/(AA + kmrgd(1))^2 , 0 ];
tens34=[ 0 , 0 ; 0 , (24*DD*kmrgd(2))/(BB^2 + kmrgd(2)^2)^2 - (72*DD*kmrgd(2)^3)/(BB^2 + kmrgd(2)^2)^3 + (48*DD*kmrgd(2)^5)/(BB^2 + kmrgd(2)^2)^4 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,2,1) =tens33;
tens3(:,:,2,2) =tens34;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,RR,AA,BB,CC,DD)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,RR,AA,BB,CC,DD)
