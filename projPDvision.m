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
f = 100e-3;
beta = 1e6;
params = [J11;J12;J13;J22;J23;J33;I22;d;m2;g;D11;D22;f;beta];

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

%Initialize more stuff.
xp(1:length(t))=0;
yp(1:length(t))=0;

%ball trajectory
% sigma(1:length(t)) = 10;
% psi = pi/4.*sin(t);
% z(1:length(t)) = 5.*sin(t);
sigma(1:length(t)) = 10;
psi(1:length(t)) = pi/400;
z(1:length(t)) = 0.01;
balltraj = [sigma; psi; z];

%control gains
Kp = 0.050.*eye(2);
Kd = 15.*eye(2);
Ki = 0.050.*eye(2);

%initialize Ei
Ei = [0; 0];

states = IC;
%RK4
for index = 1:length(t)
    q1(index) = states(1);
    q2(index) = states(2);
    q1dot(index) = states(3);
    q2dot(index) = states(4);
    
    %calculate position of ball in image
    d1 = sigma(index).*cos(q1(index)-psi(index));
    d2 = sqrt(d1.^2 + (z(index)).^2).*cos(q2(index)-atan2(z(index),sigma(index)));
    v = f.*beta./(d2 - f);
    xi = d1.*tan(q1(index)-psi(index));
    yi = d2.*tan(q2(index)-atan2(z(index),sigma(index)));
    xp(index) = v.*xi;
    yp(index) = v.*yi;
    ballimg = [xp(index); yp(index)];

    %some stuff for integral control: integral of error
    if index~=1
        Ei(1,1) = trapz(t(1:index),xp(1:index));
        Ei(2,1) = trapz(t(1:index),yp(1:index));
    end
    
    %RK4 for angles q
    k1q = derivsPDvision(t(index),states,params,ballimg,Kp,Kd,Ki,Ei);
    k2q = derivsPDvision(t(index)+h./2,states+(1/2).*h.*k1q,params,ballimg,Kp,Kd,Ki,Ei);
    k3q = derivsPDvision(t(index)+h./2,states+(1/2).*h.*k2q,params,ballimg,Kp,Kd,Ki,Ei);
    k4q = derivsPDvision(t(index)+h,states+h.*k3q,params,ballimg,Kp,Kd,Ki,Ei);
    states = states + (h./6).*(k1q + 2.*k2q + 2.*k3q + k4q);
end

figure;
hold on;
plot(t,q1,'k.-','LineWidth',3);
plot(t,q2,'b.-','LineWidth',3);
xlabel('Time (s)');
ylabel('Angles (rad)');
legend('q_{1}','q_{2}');
title('Angles Over Time: PD Control');

figure;
hold on;
plot(t,xp,'k.-','LineWidth',3);
plot(t,yp,'b.-','LineWidth',3);
legend('x_{p}','y_{p}');
xlabel('Time (s)');
ylabel('Ball Position (pix)');
title('Ball Position in Image over Time');



