%Robotics Project

clear all;
close all;

%Simulation parameters
h = 1e-3;
tend = 10;
t = 0:h:tend;

%physical parameters
J11 = 1;
J12 = 1;
J13 = 1;
J22 = 1;
J23 = 1;
J33 = J22;
I22 = 1;
d = 1;
m2 = 1;
g = 9.8;
D11 = 1;
D22 = 1;
params = [J11;J12;J13;J22;J23;J33;I22;d;m2;g;D11;D22];
D = [D11 0; 0 D22];

%initial conditions
q10 = -pi/2;
q20 = -pi/2;
q1dot0 = 0;
q2dot0 = 0;
IC = [q10;q20;q1dot0;q2dot0];

%Initialize states.
q1(1:length(t))=0;
q2(1:length(t))=0;
q1dot(1:length(t))=0;
q2dot(1:length(t))=0;

%desired trajectories
q1d = sin(t);
q2d = cos(t);
qd = [q1d; q2d];
qd1dot = cos(t);
qd2dot = -sin(t);
qd_dot = [qd1dot; qd2dot];
qd1ddot = -sin(t);
qd2ddot = -cos(t);
qd_ddot = [qd1ddot; qd2ddot];

%control gains
P = eye(9);
Kd = 10.*eye(2);
lambda = 10.*eye(2);

%Initialize ahat
ahat(1:9,1) = 0;

    %Real values for a
    areal(1,1) = J12;
    areal(2,1) = J11 + I22;
    areal(3,1) = m2.*d.^2 + J22 - J11;
    areal(4,1) = J13;
    areal(5,1) = J23;
    areal(6,1) = m2.*d.^2 + J33;
    areal(7,1) = m2.*g.*d;
    areal(8,1) = D11;
    areal(9,1) = D22;

    %Initialize some stuff for testing.
    Ya(1:2,1:length(t))=0;
    shouldbeYa(1:2,1:length(t))=0;

states = IC;
%RK4
for index = 1:length(t)
    q1(index) = states(1);
    q2(index) = states(2);
    q1dot(index) = states(3);
    q2dot(index) = states(4);

    %some initial calculations for the adaptation
    q = [q1(index); q2(index)];
    qdot = [q1dot(index); q2dot(index)];
    qtilda = q - qd(:,index);
    qtildadot = qdot - qd_dot(:,index);
    qrdot = qd_dot(:,index) - lambda*qtilda;
    qrddot = qd_ddot(:,index) - lambda*qtildadot;
    s = qdot - qrdot;
    
    %Calculate matrix Y of known stuff.
    Y(1,1) = 2.*cos(q2(index)).*sin(q2(index)).*qrddot(1) - (1 - 2.*(cos(q2(index))).^2).*(q2dot(index).*qrdot(1) + q1dot(index).*qrdot(2));
    Y(1,2) = qrddot(1);
    Y(1,3) = (cos(q2(index))).^2.*qrddot(1) - cos(q2(index)).*sin(q2(index)).*(q2dot(index).*qrdot(1) + q1dot(index).*qrdot(2));
    Y(1,4) = sin(q2(index)).*qrddot(2) + cos(q2(index)).*q2dot(index).*qrdot(2);
    Y(1,5) = cos(q2(index)).*qrddot(2) - sin(q2(index)).*q2dot(index).*qrdot(2);
    Y(1,6) = 0;
    Y(1,7) = 0;
    Y(1,8) = qrdot(1);
    Y(1,9) = 0;
    Y(2,1) = (1 - 2.*(cos(q2(index))).^2).*q1dot(index).*qrdot(1);
    Y(2,2) = 0;
    Y(2,3) = cos(q2(index)).*sin(q2(index)).*q1dot(index).*qrdot(1);
    Y(2,4) = sin(q2(index)).*qrddot(1);
    Y(2,5) = cos(q2(index)).*qrddot(1);
    Y(2,6) = qrddot(2);
    Y(2,7) = cos(q2(index));
    Y(2,8) = 0;
    Y(2,9) = qrdot(2);
    
        %for testing
        Ya(1:2,index) = Y*areal;

        H(1,1) = 2.*J12.*cos(q2(index)).*sin(q2(index)) + J11 + I22 + (m2.*d.^2 + J22 - J11).*(cos(q2(index))).^2;
        H(1,2) = J13.*sin(q2(index)) + J23.*cos(q2(index));
        H(2,1) = H(1,2);
        H(2,2) = m2.*d.^2 + J33;

        c = (1-2.*(cos(q2(index))).^2).*J12 + (m2.*d.^2 + J22 - J11).*cos(q2(index)).*sin(q2(index));
        C(1,1) = -c.*q2dot(index);
        C(1,2) = -c.*q1dot(index) + (J13.*cos(q2(index)) - J23.*sin(q2(index))).*q2dot(index);
        C(2,1) = c.*q1dot(index);
        C(2,2) = 0;

        g1 = 0;
        g2 = m2.*g.*d.*cos(q2(index));
        Gmat = [g1; g2];

        shouldbeYa(1:2,index) = H*qrddot + C*qrdot + D*qrdot + Gmat;
        %end testing block
    
    %RK4 for ahat
    k1a = ahatderivs(t(index),ahat,P,Y,s);
    k2a = ahatderivs(t(index)+h./2,ahat+(1/2).*h.*k1a,P,Y,s);
    k3a = ahatderivs(t(index)+h./2,ahat+(1/2).*h.*k2a,P,Y,s);
    k4a = ahatderivs(t(index)+h,ahat+h.*k3a,P,Y,s);
    ahat = ahat + (h./6).*(k1a + 2.*k2a + 2.*k3a + k4a);    
    
    %RK4 for angles q
    k1q = derivsAdaptive(t(index),states,params,Kd,lambda,qd(:,index),qd_dot(:,index),qd_ddot(:,index),Y,ahat);
    k2q = derivsAdaptive(t(index)+h./2,states+(1/2).*h.*k1q,params,Kd,lambda,qd(:,index),qd_dot(:,index),qd_ddot(:,index),Y,ahat);
    k3q = derivsAdaptive(t(index)+h./2,states+(1/2).*h.*k2q,params,Kd,lambda,qd(:,index),qd_dot(:,index),qd_ddot(:,index),Y,ahat);
    k4q = derivsAdaptive(t(index)+h,states+h.*k3q,params,Kd,lambda,qd(:,index),qd_dot(:,index),qd_ddot(:,index),Y,ahat);
    states = states + (h./6).*(k1q + 2.*k2q + 2.*k3q + k4q);
end

%Plot angles through time
figure;
hold on;
plot(t(1:100:end),q1d(1:100:end),'ko');
plot(t,q1,'k.-','LineWidth',3);
plot(t(1:100:end),q2d(1:100:end),'bo');
plot(t,q2,'b.-','LineWidth',3);
xlabel('Time (s)');
ylabel('Angles (rad)');
legend('q_{1d}','q_{1}','q_{2d}','q_{2}');
title('Angles Over Time: Adaptive');

%Check Y*ahat
figure;
hold on;
title('should be Ya');
plot(t,Ya(1,:),'ko','LineWidth',2);
plot(t,shouldbeYa(1,:),'k.-','LineWidth',2);
plot(t,Ya(2,:),'bo','LineWidth',2);
plot(t,shouldbeYa(2,:),'b.-','LineWidth',2);
legend('Ya row 1','shouldbeYa row 1','Ya row 2','shouldbeYa row 2');



