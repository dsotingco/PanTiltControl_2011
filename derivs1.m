%prob2qderivs.m

function res = derivs1(t,states,params)
    %extract params
    J11 = params(1);
    J12 = params(2);
    J13 = params(3);
    J22 = params(4);
    J23 = params(5);
    J33 = params(6);
    I22 = params(7);
    d = params(8);
    m2 = params(9);
    g = params(10);
    D11 = params(11);
    D22 = params(12);
    
    D = [D11 0; 0 D22];
    
    %extract states
    q1 = states(1);
    q2 = states(2);
    q = [q1; q2];
    q1dot = states(3);
    q2dot = states(4);
    qdot = [q1dot; q2dot];
    
    H(1,1) = 2.*J12.*cos(q2).*sin(q2) + J11 + I22 + (m2.*d.^2 + J22 - J11).*(cos(q2)).^2;
    H(1,2) = J13.*sin(q2) + J23.*cos(q2);
    H(2,1) = H(1,2);
    H(2,2) = m2.*d.^2 + J33;
    
    c = (1-2.*(cos(q2)).^2).*J12 + (m2.*d.^2 + J22 - J11).*cos(q2).*sin(q2);
    C(1,1) = -c.*q2dot;
    C(1,2) = -c.*q1dot + (J13.*cos(q2) - J23.*sin(q2)).*q2dot;
    C(2,1) = c.*q1dot;
    C(2,2) = 0;
    
    g1 = 0;
    g2 = m2.*g.*d.*cos(q2);
    Gmat = [g1; g2];
            
    tau = 0;
    boxterm = tau - C*qdot - D*qdot - Gmat;
    qddot = H\boxterm;
    
    res = [qdot;qddot];
end