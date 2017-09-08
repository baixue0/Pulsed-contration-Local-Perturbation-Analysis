function [splitlines] = plotStability( x,f,paramid )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
splitlines = {};
typevec = real(f(paramid,:))>0;
cross0 = find(diff(typevec));
if size(cross0,2)==0
    splitlines{1,1} = x;
    splitlines{1,2} = typevec(1)>0;
else
    cross0 = cat(2,1,cross0,size(f,2));
    for i=1:size(cross0,2)-1
        is=cross0(i); ie=cross0(i+1);
        splitlines{i,1} = x(:,is:ie);
        splitlines{i,2} = typevec(floor((is+ie)/2))>0;
    end
end

end

