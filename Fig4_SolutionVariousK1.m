% This Matlab code prints out the solutions of the system 
% k1, k2, k3 and k4
% Initialization of the parameter of the model
global k1 k2 k3 k4 C_1 C_2
k1 = .0112;
k10 =k1;
k2 = .07;
k3 = .1;
k4 = .01;
C_1 = .4;
C_2 = .8;
t0 = 0;
tf = 250;
dt = 0.01;
tspan = t0:tf;
X0 = [0.0036,0.0,1.1750,0.0]';
% Calling ode45 (RK4 method) to compute for the original (nominal)
% behaviour
[t,Y0] = ode45(@SickleCell4DE,tspan,X0);
figure(1); hold on;
subplot(2,2,1); hold on; grid on; plot(t,Y0(:,1),'b-','LineWidth',2);
subplot(2,2,2); hold on; grid on; plot(t,Y0(:,2),'b-','LineWidth',2);
subplot(2,2,3); hold on; grid on; plot(t,Y0(:,3),'b-','LineWidth',2);
subplot(2,2,4); hold on; grid on; plot(t,Y0(:,4),'b-','LineWidth',2);
k1=2*k10;
[t,Y0] = ode45(@SickleCell4DE,tspan,X0);
figure(1); 
subplot(2,2,1); plot(t,Y0(:,1),'r-','LineWidth',2);
subplot(2,2,2); plot(t,Y0(:,2),'r-','LineWidth',2);
subplot(2,2,3); plot(t,Y0(:,3),'r-','LineWidth',2);
subplot(2,2,4); plot(t,Y0(:,4),'r-','LineWidth',2);
k1=3*k10;
[t,Y0] = ode45(@SickleCell4DE,tspan,X0);
figure(1); 
subplot(2,2,1); plot(t,Y0(:,1),'g-','LineWidth',2);
subplot(2,2,2); plot(t,Y0(:,2),'g-','LineWidth',2);
subplot(2,2,3); plot(t,Y0(:,3),'g-','LineWidth',2);
subplot(2,2,4); plot(t,Y0(:,4),'g-','LineWidth',2);
k1=4*k10;
[t,Y0] = ode45(@SickleCell4DE,tspan,X0);
figure(1); 
subplot(2,2,1); plot(t,Y0(:,1),'k-','LineWidth',2);
subplot(2,2,2); plot(t,Y0(:,2),'k-','LineWidth',2);
subplot(2,2,3); plot(t,Y0(:,3),'k-','LineWidth',2);
subplot(2,2,4); plot(t,Y0(:,4),'k-','LineWidth',2);
k1=5*k10;
[t,Y0] = ode45(@SickleCell4DE,tspan,X0);
figure(1); 
subplot(2,2,1); plot(t,Y0(:,1),'m-','LineWidth',2);
subplot(2,2,2); plot(t,Y0(:,2),'m-','LineWidth',2);
subplot(2,2,3); plot(t,Y0(:,3),'m-','LineWidth',2);
subplot(2,2,4); plot(t,Y0(:,4),'m-','LineWidth',2);

figure(1); 
subplot(2,2,1); xlabel('Time [min]'); ylabel('Cm_d [mM]'); 
box on; grid on; axis([0 250 0 0.25]); %set(gca,"FontSize",10);
txt={'(a)'}; text(max(xlim)*0.1,max(ylim)*0.9,txt);
subplot(2,2,2); xlabel('Time [min]'); ylabel('Cm^{CO} [mM]'); 
box on; grid on; axis([0 250 0 1]); %set(gca,"FontSize",10);
txt={'(b)'}; text(max(xlim)*0.1,max(ylim)*0.9,txt);
subplot(2,2,3); xlabel('Time [min]'); ylabel('Cp_d [mM]'); 
box on; grid on; axis([0 250 0 1.2]);  %set(gca,"FontSize",10);
txt={'(c)'}; text(max(xlim)*0.1,max(ylim)*0.9,txt);
subplot(2,2,4); xlabel('Time [min]'); ylabel('Cp^{CO} [mM]'); 
box on; grid on; axis([0 250 0 0.4]); %set(gca,"FontSize",10);
txt={'(d)'}; text(max(xlim)*0.1,max(ylim)*0.9,txt);

ax=figure(1);
exportgraphics(ax,'Fig4.png','Resolution',1200);
