function vgdw()
clc
global mH m L a b g phi
L = 1;a = 0.5;b = L-a;
m = 5;mH = 10;g = 9.81;
options = odeset('RelTol',1e-9,'AbsTol',1e-9,'Event',@collision);
phi = 4.3*pi/180;
% Y0 = [-0.150560 0.150560 0.723611 0.586909 0]';
Y0 = [-0.306374 0.306374 1.151793 0.166594 0]';% 4.3 deg
% Y0 = [-0.298898 0.298898 1.141290 0.216675 0]';% 4 deg
% Y0 = [-0.285604 0.285604 1.120098 0.298595 0]';% 3.5 deg
% Y0 = [-0.271029 0.271029 1.093133 0.377633 0]';% 3 deg
% Y0 = [-0.254786 0.254786 1.058451 0.452139 0]';% 2.5 deg
% Y0 = [-0.236263 0.236263 1.012947 0.519347 0]';% 2 deg
% Y0 = [-0.214392 0.214392 0.951143 0.574173 0]';% 1.5 deg
% Y0 = [-0.186992 0.186992 0.861780 0.605786 0]';% 1 deg
% Y0 = [-0.148061 0.148061 0.713632 0.583338 0]';% 0.5 deg
fprintf('%f %f %f %f %f\n',Y0);
tspan = 0:0.01:2;
Y1 = [];
for i = 1:200
    [T,Y] = ode45(@eqm,tspan,Y0,options);
%     Y = [Y;Y1];
%     Y1 = Y;
    alpha = (Y(end,1)-Y(end,2))/2;
    [Qp,Qm] = tmats(alpha);
    thp = Qp\Qm*Y(end,3:4)';
    thm = Y(end,3:4)';Eloss = hsloss(thm,thp,alpha);
    Y0 = [Y(end,2) Y(end,1) thp' 0];
    fprintf('%f %f %f %f %f\n',Y0);
end
alpha=(Y0(2)-Y0(1))/2;
steplength = 2*sin(alpha)*L;
speed = steplength/T(end);
CoTmech = Eloss/(2*m+mH)/steplength/g;
fprintf('Step: %6.4f m, speed: %5.3f m/s\nEnergy Loss: %5.3f J u^2: %5.3f\n',steplength, speed, Eloss,Y(end,5)/T(end));
fprintf('CoTmech: %f \n',CoTmech);
save('vgdwdata.mat','Y')
% for i = 1:length(T)
%     [dy,u] = eqm(T(i),Y(i,:)');
%     cu(i,1:2) = u';
% end
% plot(T,cu(:,1),T,cu(:,2));
% % plot(T,Y(:,3),T,Y(:,4))
end

function [dy,u] = eqm(t,y)
global mH m L a b g phi
th1 = y(1);th2 = y(2);th1d = y(3);th2d = y(4);
M = [mH*L^2+m*L^2+m*a^2 -m*b*L*cos(th1-th2);
    -m*b*L*cos(th1-th2) m*b^2];
C = [0 -m*b*L*sin(th1-th2)*th2d;
    m*b*L*sin(th1-th2)*th1d 0];
G = [-(mH*L+m*L+m*a)*sin(th1);m*b*sin(th2)]*g;
tau = [(mH*L+m*L+m*a)*cos(th1);-m*b*cos(th2)]*g*tan(phi);
u = [1 1;0 -1]*tau;
acc = M\(tau-G-C*y(3:4,1));
dy(1:2,1) = y(3:4,1);
dy(3:4,1) = acc;
dy(5,1) = u'*u;
end

function [value,isterminal,direction] = collision(t,y)
global L
value = 1;
if y(1)> 0 && y(2) <0 && y(3)>0 && y(4)>0 && (-sin(y(1))*y(3)+sin(y(2))*y(4))<0
    value = L*(cos(y(1))-cos(y(2)));
end
isterminal = 1;
direction = -1;
end

function [Qp,Qm] = tmats(alpha)
global m mH a b L
Qp = [mH*L^2+m*a^2+m*L*(L-b*cos(2*alpha)) m*b*(b-L*cos(2*alpha));
    -m*b*L*cos(2*alpha) m*b^2];
Qm = [(mH*L^2+2*m*a*L)*cos(2*alpha)-m*a*b -m*a*b;
    -m*a*b 0];
end

function Eloss = hsloss(thdm,thdp,alpha)
global m mH L a b
M = [mH*L^2+m*L^2+m*a^2 -m*b*L*cos(2*alpha);
    -m*b*L*cos(2*alpha) m*b^2];
Eloss = (1/2)*thdm'*M*thdm-(1/2)*thdp'*M*thdp;
end
