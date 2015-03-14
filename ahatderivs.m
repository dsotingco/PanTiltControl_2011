%ahatderivs

function res = ahatderivs(t,ahat,P,Y,s)
    res = -P*Y'*s; %this is ahatdot
end

