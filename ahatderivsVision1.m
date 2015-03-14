%ahatderivs

function res = ahatderivsVision1(t,ahat,P,Y,Kp,xtilda)
    res = P*Y'*Kp*xtilda; %this is ahatdot
end

