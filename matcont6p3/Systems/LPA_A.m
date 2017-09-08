function out = LPA_A
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
function dydt = fun_eval(t,kmrgd,k0,k1,k2,k3,kT)
dydt=[(k0+k1*kmrgd(1)^3/(kmrgd(1)^3+k2^3))*(kT-kmrgd(2))-k3*kmrgd(1);;
(k0+k1*kmrgd(2)^3/(kmrgd(2)^3+k2^3))*(kT-kmrgd(2))-k3*kmrgd(2);;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(LPA_A);
y0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,k0,k1,k2,k3,kT)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,k0,k1,k2,k3,kT)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,k0,k1,k2,k3,kT)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,k0,k1,k2,k3,kT)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,k0,k1,k2,k3,kT)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,k0,k1,k2,k3,kT)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,k0,k1,k2,k3,kT)
