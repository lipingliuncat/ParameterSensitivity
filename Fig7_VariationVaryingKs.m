% This Matlab code analyzes sensitivity of the system to four parameters
% k1, k2, k3 and k4
% Initialization of the parameter of the model
global k1 k2 k3 k4 C_1 C_2
k1 = .028;
k2 = .07;
k3 = .1;
k4 = .01;
C_1 = .4;
C_2 = .1;
k0 = [k1,k2,k3,k4];
eps = min([k1,k2,k3,k4])/10.0;
t0 = 0;
tf = 400;
dt = 0.01;
tspan = t0:dt:tf;
X0 = [0.0036,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,...
    0.1750,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0]';
% Calling ode45 (RK4 method) to compute for the original (nominal)
% behaviour
[t,Y0] = ode45(@SickleCellDE,tspan,X0);

% Define the number of points for the simulation
dim = 100;
% Defines the initial starting point: 1/2 of the nominal value
k1_init = 0.5*k1;
k2_init = 0.5*k2;
k3_init = 0.5*k3;
k4_init = 0.5*k4;
% Defines the end point: 2 times of the nominal value 
k1_fin = 2.0*k1; 
k2_fin = 2.0*k2; 
k3_fin = 2.0*k3; 
k4_fin = 2.0*k4; 
% Creates a uniform random number distribution between the set of 2 points
k1v = random('unif',k1_init,k1_fin,1,dim);
k2v = random('unif',k2_init,k2_fin,1,dim);
k3v = random('unif',k3_init,k3_fin,1,dim);
k4v = random('unif',k4_init,k4_fin,1,dim);
k1 = k1v(1);
k2 = k2v(1);
k3 = k3v(1);
k4 = k4v(1);
[t,Y] = ode45(@SickleCellDE,tspan,X0);
figure(1); 
subplot(2,2,1); hold on; grid on; plot(t,Y(:,1)-Y0(:,1),'LineWidth',1);
subplot(2,2,2); hold on; grid on; plot(t,Y(:,6)-Y0(:,6),'LineWidth',1);
subplot(2,2,3); hold on; grid on; plot(t,Y(:,11)-Y0(:,11),'LineWidth',1);
subplot(2,2,4); hold on; grid on; plot(t,Y(:,16)-Y0(:,16),'LineWidth',1);

for i_k = 2:1:dim
    k1 = k1v(i_k);
    k2 = k2v(i_k);
    k3 = k3v(i_k);
    k4 = k4v(i_k);
    [t,Y] = ode45(@SickleCellDE,tspan,X0);
    
    figure(1); 
    subplot(2,2,1);     plot(t,Y(:,1)-Y0(:,1),'LineWidth',1);
    subplot(2,2,2);     plot(t,Y(:,6)-Y0(:,6),'LineWidth',1);
    subplot(2,2,3);     plot(t,Y(:,11)-Y0(:,11),'LineWidth',1);
    subplot(2,2,4);     plot(t,Y(:,16)-Y0(:,16),'LineWidth',1);
     
end
% Add titles to figures in for loop1 above
figure(1); 
subplot(2,2,1); 
box on; grid on; xlabel('Time [min]'); 
ylabel('$x(t)-x(t)_\mathrm{ref}$','Interpreter','latex');
txt={'(a)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);
subplot(2,2,2);  
box on; grid on; xlabel('Time [min]'); 
ylabel('$y(t)-y(t)_\mathrm{ref}$','Interpreter','latex');
txt={'(b)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);
subplot(2,2,3); 
box on; grid on; xlabel('Time [min]'); 
ylabel('$z(t)-z(t)_\mathrm{ref}$','Interpreter','latex');
txt={'(c)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);
subplot(2,2,4); 
box on; grid on; xlabel('Time [min]'); 
ylabel('$u(t)-u(t)_\mathrm{ref}$','Interpreter','latex');
txt={'(d)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);

ax=figure(1);
exportgraphics(ax,'Fig7.png','Resolution',1200);
