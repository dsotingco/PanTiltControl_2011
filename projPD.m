%Robotics Project

clear all;
close all;

%Simulation parameters
h = 5e-2;
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

%initial conditions
q10 = 0;
q20 = 0;
q1dot0 = 0;
q2dot0 = 0;
IC = [q10;q20;q1dot0;q2dot0];

%Initialize states.
q1(1:length(t))=0;
q2(1:length(t))=0;
q1dot(1:length(t))=0;
q2dot(1:length(t))=0;

%desired angles
q1d = pi/2;
q2d = pi;
qd = [q1d; q2d];

%control gains
Kp = 100.*eye(2);
Kd = 10.*eye(2);

states = IC;
%RK4
for index = 1:length(t)
    q1(index) = states(1);
    q2(index) = states(2);
    q1dot(index) = states(3);
    q2dot(index) = states(4);

    %RK4 for angles q
    k1q = derivsPD(t(index),states,params,qd,Kp,Kd);
    k2q = derivsPD(t(index)+h./2,states+(1/2).*h.*k1q,params,qd,Kp,Kd);
    k3q = derivsPD(t(index)+h./2,states+(1/2).*h.*k2q,params,qd,Kp,Kd);
    k4q = derivsPD(t(index)+h,states+h.*k3q,params,qd,Kp,Kd);
    states = states + (h./6).*(k1q + 2.*k2q + 2.*k3q + k4q);
end

figure;
hold on;
plot([0 tend],[q1d q1d],'k.--','LineWidth',2);
plot(t,q1,'k.-','LineWidth',3);
plot([0 tend],[q2d q2d],'b.--','LineWidth',2);
plot(t,q2,'b.-','LineWidth',3);
xlabel('Time (s)');
ylabel('Angles (rad)');
legend('q_{1d}','q_{1}','q_{2d}','q_{2}');
title('Angles Over Time: PD Control');



