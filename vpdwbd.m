function vpdwbd()
clc
global mH m L a b g phi
L = 1;a = 0.5;b = L-a;
m = 5;mH = 10;g = 9.81;
options = odeset('RelTol',1e-9,'AbsTol',1e-9,'Event',@collision);
Y0 = [-0.050358 0.050358 0.258993 0.253712]';
sol = zeros(200,9);
Yprev = zeros(4,1);
for j = 1:1:225
    phi = 2*j/100*pi/180;
    sol(j,1) = phi;
%     fprintf('%f %f %f %f\n',Y0);
    tic
    tspan = 0:0.01:2;
    for i = 1:100
        [T,Y] = ode45(@eqm,tspan,Y0,options);
        alpha = (Y(end,1)-Y(end,2))/2;
        [Qp,Qm] = tmats(alpha);
        thp = Qp\Qm*Y(end,3:4)';
        thm = Y(end,3:4)';Eloss = hsloss(thm,thp,alpha);
        Yprev = Y0;
        Y0 = [Y(end,2) Y(end,1) thp']';
    %     fprintf('%f %f %f %f\n',Y0);
    end
 
    fprintf('Phi: %f deg, %f %f %f %f\n',phi*180/pi,Y0);
    toc
    sol(j,2:9) = [Y0' Yprev'];
end
save('vpdwbddata.mat','sol')
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
acc = M\(tau-G-C*y(3:4,1));
dy(1:2,1) = y(3:4,1);
dy(3:4,1) = acc;
end

function [value,isterminal,direction] = collision(t,y)
global L
value = 1;
if y(1)> 0 && y(2) <0 && y(3)>0 && y(4)>0
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
