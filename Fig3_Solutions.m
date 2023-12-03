% This Matlab code prints out the solutions of the system 
% k1, k2, k3 and k4
% Initialization of the parameter of the model
global k1 k2 k3 k4 C_1 C_2
k1 = .0112;
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
figure(1); hold on; grid on;
subplot(2,2,1); hold on; grid on; plot(t,Y0(:,1),'b-','LineWidth',2);
subplot(2,2,2); hold on; grid on; plot(t,Y0(:,2),'r-','LineWidth',2);
plot(t,Y0(:,3),'g-','LineWidth',2);
plot(t,Y0(:,4),'k-','LineWidth',2);
%legend('De-oxy Monomers [x]','CO-bound Monomers [y]',...
%    'De-oxy polymers [z]','CO-bound Polymers [u]' );
xlabel('Time [min]'); %  label the x-axis
ylabel('Molar Concentration [mM]'); %  label the y-axis
box on;
set(gca,"FontSize",12);
ax=gca;
exportgraphics(ax,'Fig3PNGr.png','Resolution',1200);
