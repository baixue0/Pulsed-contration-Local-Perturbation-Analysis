function out = m2
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9)
dydt=[(k0+k1*kmrgd(1)^3/(kmrgd(1)^3+k2^3))*(kmrgd(3)-kmrgd(2))-(k3+k4*kmrgd(4)/(1+kmrgd(4)))*kmrgd(1);
(k0+k1*kmrgd(2)^3/(kmrgd(2)^3+k2^3))*(kmrgd(3)-kmrgd(2))-(k3+k4*kmrgd(5)/(1+kmrgd(5)))*kmrgd(2);
k8-k9*kmrgd(3);%total A+I on membrane
k5*kmrgd(1)-k6*kmrgd(4);
k5*kmrgd(2)-k6*kmrgd(5);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(m2);
y0=[0,0,0,0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9)
