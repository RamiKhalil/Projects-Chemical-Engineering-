function fp=der(a,x,h)
% Calculates the derivate of a function numerically based
% on the 5-pt Stencil
% (c) 2011 Phillip Servio
fp=(-feval(a,x+2*h)+8*feval(a,x+h)-8*feval(a,x-h)...
    +feval(a,x-2*h))/(12*h);
end
