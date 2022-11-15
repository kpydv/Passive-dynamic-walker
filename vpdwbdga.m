function vpdwbdga()
global mH m L a b g phi
L = 1;a = 0.5;b = L-a;
m = 5;mH = 10;g = 9.81;phi = 2*pi/180;
clc
options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', ...
    @gaplotscores,'MaxGenerations',10000);
[x,fval,exitflag,output,population,scores] = ga(@vpdwfitness,3,[],[],[],[],[0.1 0.1 0.1],[0.5 2 2],[],options)
% [X,FVAL] = fmincon(@vpdwfitness,[0.138176 0.583781 0.174595],[],[],[],[],[0.1 0 0],[1 1 1])
format long
[-x(1) x]'
format

end

function err = vpdwfitness(par)
Y0 = [-par(1) par]';
options = odeset('RelTol',1e-4,'AbsTol',1e-4,'Event',@collision);
tspan = 0:0.01:2;
[T,Y] = ode45(@eqm,tspan,Y0,options);
alpha = Y(end,1)-Y(end,2);
[Qp,Qm] = tmats(alpha);
thp = Qp\Qm*Y(end,3:4)';
Y0new = [Y(end,2) Y(end,1) thp'];
err = norm(Y0new'-Y0);
end

function [dy] = eqm(t,y)
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
if y(1)>0 || y(1)-y(2)>0.1
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

function compassAnim(t,x)
global l phi
gam = phi;
th_st = x(1,1); th_sw = x(1,2);
% postion of stance leg
xst = -1.051;
yst = 0;
% position of hip wrt to stance point
xhp = xst - l*sin(th_st(1));
yhp = yst + l*cos(th_st(1));
% position of swing leg foot point
xsw = l*sin( th_sw(1))+ xhp;
ysw = yhp - l*cos( th_sw(1));
axis([-1+xhp 1+xhp -1 1.5])
slope = line([-2.23 10.5],[0.042+0.09 (-2.23-13.5)*tan(gam)],'Color','k','LineWidth',1);
set(slope,'xdata',[-2.23 10.5],'ydata',[0.042+0.09 (-2.23-13.5)*tan(gam)]);
stleg = line([xst(1) xhp(1)],[yst(1) yhp(1)],'Color','k','LineWidth',2);
swleg = line([xsw(1) xhp(1)],[ysw(1) yhp(1)],'Color','b','LineWidth',2);

for j = 1:length(x(:,1))
th_st = x(j,1); th_sw = x(j,2);
xst = -.51;
yst = 0;
xhp = xst - l*sin(th_st);
yhp = yst + l*cos(th_st);
xsw = l*sin( th_sw)+ xhp;
ysw = yhp - l*cos(th_sw);
axis([-1+xhp 1+xhp -1 1.5])
set(stleg,'xdata',[xst xhp],'ydata',[yst yhp]);
set(swleg,'xdata',[xsw xhp],'ydata',[ysw yhp]);

drawnow
pause(0.02);
axis off
hold off
end
end