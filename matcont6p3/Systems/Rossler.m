function out = Rossler
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,AA,BB,CC)
dydt=[-kmrgd(2)-kmrgd(3);;
kmrgd(1)+AA*kmrgd(2);;
BB*kmrgd(1)-CC*kmrgd(3)+kmrgd(1)*kmrgd(3);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(Rossler);
y0=[0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,AA,BB,CC)
jac=[ 0 , -1 , -1 ; 1 , AA , 0 ; BB + kmrgd(3) , 0 , kmrgd(1) - CC ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,AA,BB,CC)
jacp=[ 0 , 0 , 0 ; kmrgd(2) , 0 , 0 ; 0 , kmrgd(1) , -kmrgd(3) ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,AA,BB,CC)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,AA,BB,CC)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,AA,BB,CC)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,AA,BB,CC)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,AA,BB,CC)
