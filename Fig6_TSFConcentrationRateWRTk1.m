global k1 k2 k3 k4 C_1 C_2
k1 = .028;
k2 = .07;
k3 = .1;
k4 = .01;
C_1 = .4;
C_2 = .1;
k0 = [k1,k2,k3,k4];
t0 = 0;
tf = 400;
tspan = t0:tf;
X0 = [0.0036,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,...
    0.1750,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0]';
[t,Y0] = ode45(@SickleCellDE,tspan,X0);
figure(1); 
subplot(2,2,1); hold on; grid on; plot(t,Y0(:,2),'b-','LineWidth',2);
box on; grid on; xlabel('Time [min]'); 
ylabel('$\partial x(t) / \partial k_1$','Interpreter','latex');
txt={'(a)'}; text(max(xlim)*0.1,max(ylim)*0.9,txt);
%set(gca,"FontSize",10);
subplot(2,2,2); hold on; grid on; plot(t,Y0(:,7),'b-','LineWidth',2);
box on; grid on; xlabel('Time [min]'); 
ylabel('$\partial y(t) / \partial k_1$','Interpreter','latex');
txt={'(b)'}; text(max(xlim)*0.1,max(ylim)*0.9,txt);
subplot(2,2,3); hold on; grid on; plot(t,Y0(:,12),'b-','LineWidth',2);
box on; grid on; xlabel('Time [min]'); 
ylabel('$\partial z(t) / \partial k_1$','Interpreter','latex');
txt={'(c)'}; text(max(xlim)*0.1,min(ylim)*0.1,txt);
subplot(2,2,4); hold on; grid on; plot(t,Y0(:,17),'b-','LineWidth',2);
box on; grid on; xlabel('Time [min]'); 
ylabel('$\partial u(t) / \partial k_1$','Interpreter','latex');
txt={'(d)'}; text(max(xlim)*0.1,max(ylim)*0.9,txt);

ax=figure(1);
exportgraphics(ax,'Fig6.png','Resolution',1200);
