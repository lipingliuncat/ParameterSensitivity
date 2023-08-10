%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this code, we are conducting the simulation for the binding and  
% melting parameters of the sickle cell model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code analyzes sensitivity of the system to four paramters
% k1, k2, k3 and k4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
X0 = [0.0036,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,...
    0.1750,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0];

% Calling My runge-kutta 4 method to compute the original (nominal)
% behaviour
[t,Y0] = MYRK4COSens(X0,t0,tf,dt);
product = Y0(:,6);
figure(1)
hold on; grid on;
plot(t,Y0(:,1),'b')
plot(t,Y0(:,6),'r')
plot(t,Y0(:,11),'g')
plot(t,Y0(:,16),'k')
legend('De-oxy Monomers [x]','CO-bound Monomers [y]',...
    'De-oxy polymers [z]','CO-bound Polymers [u]' );
title('Molar Cencentration') % put a title on the plot
xlabel('Time [min]')%  label the x-axis
ylabel('Population [Mm]') %  label the y-axis
h1 = figure(1);
saveas(h1,'Fig1_solnCurves4model','jpg')

figure(2)
hold on; grid on;
subplot(2,2,1)
hold on; grid on;
plot(t,Y0(:,2))
title('dx/dk1 solution curve');
subplot(2,2,2)
hold on; grid on;
plot(t,Y0(:,7))
title('dy/dk1 solution curve');
subplot(2,2,3)
hold on; grid on;
plot(t,Y0(:,12))
title('dz/dk1 solution curve');
subplot(2,2,4)
hold on; grid on;
plot(t,Y0(:,17))
title('du/dk1 solution curve');
suptitle('Behavior of dx/dk1,dy/dk1,dz/dk1,du/dk1 variables');
h2 = figure(2);
saveas(h2,'Fig2_solnCurves4variationsOfK1','jpg')

figure(3)
hold on; grid on;
subplot(2,2,1)
hold on; grid on;
plot(t,Y0(:,3))
title('dx/dk2 solution curve');
subplot(2,2,2)
hold on; grid on;
plot(t,Y0(:,8))
title('dy/dk2 solution curve');
subplot(2,2,3)
hold on; grid on;
plot(t,Y0(:,13))
title('dz/dk2 solution curve');
subplot(2,2,4)
hold on; grid on;
plot(t,Y0(:,18))
title('du/dk2 solution curve');
suptitle('Behavior of dx/dk2,dy/dk2,dz/dk2,du/dk2 variables');
h3 = figure(3);
saveas(h3,'Fig3_solnCurves4variationsOfK2','jpg')

figure(4)
hold on; grid on;
subplot(2,2,1)
hold on; grid on;
plot(t,Y0(:,4))
title('dx/dk3 solution curve');
subplot(2,2,2)
hold on; grid on;
plot(t,Y0(:,9))
title('dy/dk3 solution curve');
subplot(2,2,3)
hold on; grid on;
plot(t,Y0(:,14))
title('dz/dk3 solution curve');
subplot(2,2,4)
hold on; grid on;
plot(t,Y0(:,19))
title('du/dk3 solution curve');
suptitle('Behavior of dx/dk3,dy/dk3,dz/dk3,du/dk3 variables');
h4 = figure(4);
saveas(h4,'Fig4_solnCurves4variationsOfK3','jpg')

figure(5)
hold on; grid on;
subplot(2,2,1)
hold on; grid on;
plot(t,Y0(:,5))
title('dx/dk4 solution curve');
subplot(2,2,2)
hold on; grid on;
plot(t,Y0(:,10))
title('dy/dk4 solution curve');
subplot(2,2,3)
hold on; grid on;
plot(t,Y0(:,15))
title('dz/dk4 solution curve');
subplot(2,2,4)
hold on; grid on;
plot(t,Y0(:,20))
title('du/dk4 solution curve');
suptitle('Behavior of dx/dk4,dy/dk4,dz/dk4,du/dk4 variables');
h5 = figure(5);
saveas(h5,'Fig5_solnCurves4variationsOfK4','jpg')

saveas(h1,'Fig1_solnCurves4model','png')
saveas(h2,'Fig2_solnCurves4variationsOfK1','png')
saveas(h3,'Fig3_solnCurves4variationsOfK2','png')
saveas(h4,'Fig4_solnCurves4variationsOfK3','png')
saveas(h5,'Fig5_solnCurves4variationsOfK4','png')

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

[t,Y] = MYRK4COSens(X0,t0,tf,dt);

figure(6)
hold on; grid on;
subplot(2,2,1)
hold on; grid on;
plot(t,Y(:,1)-Y0(:,1))
ylabel('x-x_ref')
%title('x variation wrt to Ks');
subplot(2,2,2)
hold on; grid on;
plot(t,Y(:,6)-Y0(:,6))
ylabel('y-y_ref')
%title('y variation wrt Ks');
subplot(2,2,3)
hold on; grid on;
plot(t,Y(:,11)-Y0(:,11))
ylabel('z-z_ref')
%title('z variation wrt Ks');
subplot(2,2,4)
hold on; grid on;
plot(t,Y(:,16)-Y0(:,16))
ylabel('u-u_ref')

k1x_min = abs(k1-k0(1))*sum(abs(Y(:,2)-Y0(:,2)).^2)*dt;
k1y_min = abs(k1-k0(1))*sum(abs(Y(:,7)-Y0(:,7)).^2)*dt;
k1z_min = abs(k1-k0(1))*sum(abs(Y(:,12)-Y0(:,12)).^2)*dt;
k1u_min = abs(k1-k0(1))*sum(abs(Y(:,17)-Y0(:,17)).^2)*dt;
    
k2x_min = abs(k2-k0(2))*sum(abs(Y(:,3)-Y0(:,3)).^2)*dt;
k2y_min = abs(k2-k0(2))*sum(abs(Y(:,8)-Y0(:,8)).^2)*dt;
k2z_min = abs(k2-k0(2))*sum(abs(Y(:,13)-Y0(:,13).^2))*dt;
k2u_min = abs(k2-k0(2))*sum(abs(Y(:,18)-Y0(:,18)).^2)*dt;
    
k3x_min = abs(k3-k0(3))*sum(abs(Y(:,4)-Y0(:,4)).^2)*dt;
k3y_min = abs(k3-k0(3))*sum(abs(Y(:,9)-Y0(:,9)).^2)*dt;
k3z_min = abs(k3-k0(3))*sum(abs(Y(:,14)-Y0(:,14)).^2)*dt;
k3u_min = abs(k3-k0(3))*sum(abs(Y(:,19)-Y0(:,19)).^2)*dt;
    
k4x_min = abs(k4-k0(4))*sum(abs(Y(:,5)-Y0(:,5)).^2)*dt;
k4y_min = abs(k4-k0(4))*sum(abs(Y(:,10)-Y0(:,10)).^2)*dt;
k4z_min = abs(k4-k0(4))*sum(abs(Y(:,15)-Y0(:,15)).^2)*dt;
k4u_min = abs(k4-k0(4))*sum(abs(Y(:,20)-Y0(:,20)).^2)*dt;

k1x_max = abs(k1-k0(1))*sum(abs(Y(:,2)-Y0(:,2)).^2)*dt;
k1y_max = abs(k1-k0(1))*sum(abs(Y(:,7)-Y0(:,7)).^2)*dt;
k1z_max = abs(k1-k0(1))*sum(abs(Y(:,12)-Y0(:,12)).^2)*dt;
k1u_max = abs(k1-k0(1))*sum(abs(Y(:,17)-Y0(:,17)).^2)*dt;
    
k2x_max = abs(k2-k0(2))*sum(abs(Y(:,3)-Y0(:,3)).^2)*dt;
k2y_max = abs(k2-k0(2))*sum(abs(Y(:,8)-Y0(:,8)).^2)*dt;
k2z_max = abs(k2-k0(2))*sum(abs(Y(:,13)-Y0(:,13)).^2)*dt;
k2u_max = abs(k2-k0(2))*sum(abs(Y(:,18)-Y0(:,18)).^2)*dt;
    
k3x_max = abs(k3-k0(3))*sum(abs(Y(:,4)-Y0(:,4)).^2)*dt;
k3y_max = abs(k3-k0(3))*sum(abs(Y(:,9)-Y0(:,9)).^2)*dt;
k3z_max = abs(k3-k0(3))*sum(abs(Y(:,14)-Y0(:,14)).^2)*dt;
k3u_max = abs(k3-k0(3))*sum(abs(Y(:,19)-Y0(:,19)).^2)*dt;
    
k4x_max = abs(k4-k0(4))*sum(abs(Y(:,5)-Y0(:,5)).^2)*dt;
k4y_max = abs(k4-k0(4))*sum(abs(Y(:,10)-Y0(:,10)).^2)*dt;
k4z_max = abs(k4-k0(4))*sum(abs(Y(:,15)-Y0(:,15)).^2)*dt;
k4u_max = abs(k4-k0(4))*sum(abs(Y(:,20)-Y0(:,20)).^2)*dt;

x_diff_min = sum(abs(Y(:,1)-Y0(:,1)).^2)*dt;
y_diff_min = sum(abs(Y(:,6)-Y0(:,6)).^2)*dt;
z_diff_min = sum(abs(Y(:,11)-Y0(:,11)).^2)*dt;
u_diff_min = sum(abs(Y(:,16)-Y0(:,16)).^2)*dt;

x_diff_max = sum(abs(Y(:,1)-Y0(:,1)).^2)*dt;
y_diff_max = sum(abs(Y(:,6)-Y0(:,6)).^2)*dt;
z_diff_max = sum(abs(Y(:,11)-Y0(:,11)).^2)*dt;
u_diff_max = sum(abs(Y(:,16)-Y0(:,16)).^2)*dt;
    

% For loop1: sensitivity of x,y,z,u wrt k1,k2,k3,k4 when all four 
% parameters vary simultaneously
for i_k = 2:1:dim
    k1 = k1v(i_k);
    k2 = k2v(i_k);
    k3 = k3v(i_k);
    k4 = k4v(i_k);
    
    [t,Y] = MYRK4COSens(X0,t0,tf,dt);
    x_diff = sum(abs(Y(:,1)-Y0(:,1)).^2)*dt;
    y_diff = sum(abs(Y(:,6)-Y0(:,6)).^2)*dt;
    z_diff = sum(abs(Y(:,11)-Y0(:,11).^2))*dt;
    u_diff = sum(abs(Y(:,16)-Y0(:,16)).^2)*dt;
    
    if (x_diff<x_diff_min)
        x_diff_min = x_diff;
    elseif (x_diff>x_diff_max)
        x_diff_max = x_diff;
    end
    
    if (y_diff<y_diff_min)
        y_diff_min = y_diff;
    elseif (y_diff>y_diff_max)
        y_diff_max = y_diff;
    end
    
    if (z_diff<z_diff_min)
        z_diff_min = z_diff;
    elseif (x_diff>x_diff_max)
        z_diff_max = z_diff;
    end
    
    if (u_diff<u_diff_min)
        u_diff_min = u_diff;
    elseif (u_diff>u_diff_max)
        u_diff_max = u_diff;
    end
    
    figure(6)
    hold on; grid on;
    subplot(2,2,1)
    hold on; grid on;
    plot(t,Y(:,1)-Y0(:,1))
    ylabel('x-x_ref')
    %title('x variation wrt to Ks');
    subplot(2,2,2)
    hold on; grid on;
    plot(t,Y(:,6)-Y0(:,6))
    ylabel('y-y_ref')
    %title('y variation wrt Ks');
    subplot(2,2,3)
    hold on; grid on;
    plot(t,Y(:,11)-Y0(:,11))
    ylabel('z-z_ref')
    %title('z variation wrt Ks');
    subplot(2,2,4)
    hold on; grid on;
    plot(t,Y(:,16)-Y0(:,16))
    ylabel('u-u_ref')
    %title('u variation wrt k1');
    %suptitle('Behavior of x,y,z,u variables wrt Ks');

    k1x_tmp = abs(k1-k0(1))*sum(abs(Y(:,2)-Y0(:,2)).^2)*dt;
    k1y_tmp = abs(k1-k0(1))*sum(abs(Y(:,7)-Y0(:,7)).^2)*dt;
    k1z_tmp = abs(k1-k0(1))*sum(abs(Y(:,12)-Y0(:,12)).^2)*dt;
    k1u_tmp = abs(k1-k0(1))*sum(abs(Y(:,17)-Y0(:,17)).^2)*dt;
    
    k2x_tmp = abs(k2-k0(2))*sum(abs(Y(:,3)-Y0(:,3)).^2)*dt;
    k2y_tmp = abs(k2-k0(2))*sum(abs(Y(:,8)-Y0(:,8)).^2)*dt;
    k2z_tmp = abs(k2-k0(2))*sum(abs(Y(:,13)-Y0(:,13)).^2)*dt;
    k2u_tmp = abs(k2-k0(2))*sum(abs(Y(:,18)-Y0(:,18)).^2)*dt;
    
    k3x_tmp = abs(k3-k0(3))*sum(abs(Y(:,4)-Y0(:,4)).^2)*dt;
    k3y_tmp = abs(k3-k0(3))*sum(abs(Y(:,9)-Y0(:,9)).^2)*dt;
    k3z_tmp = abs(k3-k0(3))*sum(abs(Y(:,14)-Y0(:,14)).^2)*dt;
    k3u_tmp = abs(k3-k0(3))*sum(abs(Y(:,19)-Y0(:,19)).^2)*dt;
    
    k4x_tmp = abs(k4-k0(4))*sum(abs(Y(:,5)-Y0(:,5)).^2)*dt;
    k4y_tmp = abs(k4-k0(4))*sum(abs(Y(:,10)-Y0(:,10)).^2)*dt;
    k4z_tmp = abs(k4-k0(4))*sum(abs(Y(:,15)-Y0(:,15)).^2)*dt;
    k4u_tmp = abs(k4-k0(4))*sum(abs(Y(:,20)-Y0(:,20)).^2)*dt;

    
    if (k1x_tmp < k1x_min)
        k1x_min = k1x_tmp;
    elseif (k1x_tmp > k1x_max)
        k1x_max = k1x_tmp;
    end
    
    if (k1y_tmp < k1y_min)
        k1y_min = k1y_tmp;
    elseif (k1y_tmp > k1y_max)
        k1y_max = k1y_tmp;
    end
    
    if (k1z_tmp < k1z_min)
        k1z_min = k1z_tmp;
    elseif (k1z_tmp > k1z_max)
        k1z_max = k1z_tmp;
    end
    
    if (k1u_tmp < k1u_min)
        k1u_min = k1u_tmp;
    elseif (k1u_tmp > k1u_max)
        k1u_max = k1u_tmp;
    end
    
    if (k2x_tmp < k2x_min)
        k2x_min = k2x_tmp;
    elseif (k2x_tmp > k2x_max)
        k2x_max = k2x_tmp;
    end
    
    if (k2y_tmp < k2y_min)
        k2y_min = k2y_tmp;
    elseif (k2y_tmp > k2y_max)
        k2y_max = k2y_tmp;
    end
    
    if (k2z_tmp < k2z_min)
        k2z_min = k2z_tmp;
    elseif(k2z_tmp > k2z_max)
        k2z_max = k2z_tmp;
    end
    
    if (k2u_tmp < k2u_min)
        k2u_min = k2u_tmp;
    elseif (k2u_tmp > k2u_max)
        k2u_max = k2u_tmp;
    end
    
    if (k3x_tmp < k3x_min)
        k3x_min = k3x_tmp;
    elseif (k3x_tmp > k3x_max)
        k3x_max = k3x_tmp;
    end
    
    if (k3y_tmp < k3y_min)
        k3y_min = k3y_tmp;
    elseif (k3y_tmp > k3y_max)
        k3y_max = k3y_tmp;
    end
    
    if (k3z_tmp < k3z_min)
        k3z_min = k3z_tmp;
    elseif(k3z_tmp > k3z_max)
        k3z_max = k3z_tmp;
    end
    
    if (k3u_tmp < k3u_min)
        k3u_min = k3u_tmp;
    elseif (k3u_tmp > k3u_max)
        k3u_max = k3u_tmp;
    end
    
    if (k4x_tmp < k4x_min)
        k4x_min = k4x_tmp;
    elseif (k4x_tmp > k4x_max)
        k4x_max = k4x_tmp;
    end
    
    if (k4y_tmp < k4y_min)
        k4y_min = k4y_tmp;
    elseif (k4y_tmp > k4y_max)
        k4y_max = k4y_tmp;
    end
    
    if (k4z_tmp < k4z_min)
        k4z_min = k4z_tmp;
    elseif(k4z_tmp > k4z_max)
        k4z_max = k4z_tmp;
    end
    
    if (k4u_tmp < k4u_min)
        k4u_min = k4u_tmp;
    elseif (k4u_tmp > k4u_max)
        k4u_max = k4u_tmp;
    end
    
end

% Add titles to figures in for loop1 above
figure(6)
hold on; grid on;
subplot(2,2,1)
hold on; grid on;
ylabel('x-x_ref')
title('x variation wrt varying Ks');
subplot(2,2,2)
hold on; grid on;
ylabel('y-y_ref')
title('y variation wrt varying Ks');
subplot(2,2,3)
hold on; grid on;
ylabel('z-z_ref')
title('z variation wrt varying Ks');
subplot(2,2,4)
hold on; grid on;
ylabel('u-u_ref')
title('u variation wrt varying Ks');
suptitle('Behavior of x,y,z,u variables for varying Ks');

h1 = figure(6);
saveas(h1,'Fig6_DiffWRefSolnCurve_varyKs','jpg')
saveas(h1,'Fig6_DiffWRefSolnCurve_varyKs','png')


RefSolnDiff_range = [x_diff_min,y_diff_min,z_diff_min,u_diff_min;
    x_diff_max,y_diff_max,z_diff_max,u_diff_max;
    x_diff_max-x_diff_min,y_diff_max-y_diff_min,z_diff_max-z_diff_min,...
    u_diff_max-u_diff_min]

fileID = fopen('SolnDiffRange.txt','w');
fprintf(fileID,'%10s %12s %12s %12s\r\n','x_diff','y_diff',...
    'z_diff','u_diff');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\r\n',RefSolnDiff_range');
fclose(fileID);

k1_range = [k1x_min,k1y_min,k1z_min,k1u_min;
    k1x_max,k1y_max,k1z_max,k1u_max;
    k1x_max-k1x_min,k1y_max-k1y_min,k1z_max-k1z_min,k1u_max-k1u_min]

fileID = fopen('K1Range.txt','w');
fprintf(fileID,'%10s %12s %12s %12s\r\n','k1x','k1y','k1z','k1u');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\r\n',k1_range');
fclose(fileID);

k2_range = [k2x_min,k2y_min,k2z_min,k2u_min;
    k2x_max,k2y_max,k2z_max,k2u_max;
    k2x_max-k2x_min,k2y_max-k2y_min,k2z_max-k2z_min,k2u_max-k2u_min]

fileID = fopen('K2Range.txt','w');
fprintf(fileID,'%10s %12s %12s %12s\r\n','k2x','k2y','k2z','k2u');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\r\n',k2_range');
fclose(fileID);

k3_range = [k3x_min,k3y_min,k3z_min,k3u_min;
    k3x_max,k3y_max,k3z_max,k3u_max;
    k3x_max-k3x_min,k3y_max-k3y_min,k3z_max-k3z_min,k3u_max-k3u_min]

fileID = fopen('K3Range.txt','w');
fprintf(fileID,'%10s %12s %12s %12s\r\n','k3x','k3y','k3z','k3u');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\r\n',k3_range');
fclose(fileID);

k4_range = [k4x_min,k4y_min,k4z_min,k4u_min;
    k4x_max,k4y_max,k4z_max,k4u_max;
    k4x_max-k4x_min,k4y_max-k4y_min,k4z_max-k4z_min,k4u_max-k4u_min]

fileID = fopen('K4Range.txt','w');
fprintf(fileID,'%10s %12s %12s %12s\r\n','k4x','k4y','k4z','k4u');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\r\n',k4_range');
fclose(fileID);

% For loop1: sensitivity of x,y,z,u wrt k1,k2,k3,k4 when all four 
% parametersvary simultaneously
    for i_k = 1:1:dim
        k1 = k1v(i_k);
        k2 = k2v(i_k);
        k3 = k3v(i_k);
        k4 = k4v(i_k);

        [t,Y] = MYRK4COSensNew(X0,t0,tf,dt);


        figure(6)
        hold on; grid on;
        subplot(2,2,1)
        hold on; grid on;
        plot(t,Y(:,1))
        ylabel('x')
        title('x variation wrt to k1');
        subplot(2,2,2)
        hold on; grid on;
        plot(t,Y(:,6))
        ylabel('y')
        title('y variation wrt k1');
        subplot(2,2,3)
        hold on; grid on;
        plot(t,Y(:,11))
        ylabel('z')
        title('z variation wrt k1');
        subplot(2,2,4)
        hold on; grid on;
        plot(t,Y(:,16))
        ylabel('u')
        title('u variation wrt k1');
        suptitle('Behavior of x,y,z,u variables wrt k1');
    
        figure(7)
        hold on; grid on;
        subplot(2,2,1)
        hold on; grid on;
        plot(t,Y(:,2))
        ylabel('dx/dk1')
        title('x variation wrt to k1');
        subplot(2,2,2)
        hold on; grid on;
        plot(t,Y(:,7))
        ylabel('dy/dk1')
        title('y variation wrt k1');
        subplot(2,2,3)
        hold on; grid on;
        plot(t,Y(:,12))
        ylabel('dz/dk1')
        title('z variation wrt k1');
        subplot(2,2,4)
        hold on; grid on;
        plot(t,Y(:,17))
        ylabel('du/dk1')
        title('u variation wrt k1');
        suptitle('Behavior of dx,dy,dz,du variables wrt k1');
    
        figure(8)
        hold on; grid on;
        subplot(2,2,1)
        hold on; grid on;
        plot(t,Y(:,3))
        ylabel('dx/dk2')
        title('x variation wrt k2');
        subplot(2,2,2)
        hold on; grid on;
        plot(t,Y(:,8))
        ylabel('dy/dk2')
        title('y variation wrt k2');
        subplot(2,2,3)
        hold on; grid on;
        plot(t,Y(:,13))
        ylabel('dz/dk2')
        title('z variation wrt k2');
        subplot(2,2,4)
        hold on; grid on;
        plot(t,Y(:,18))
        ylabel('du/dk2')
        title('u variation wrt k2');
        suptitle('Behavior of x,y,z,u variables wrt k2');
    
        figure(9)
        hold on; grid on;
        subplot(2,2,1)
        hold on; grid on;
        plot(t,Y(:,4))
        ylabel('dx/dk3')
        title('x variation wrt k3');
        subplot(2,2,2)
        hold on; grid on;
        plot(t,Y(:,9))
        ylabel('dy/dk3')
        title('y variation wrt k3');
        subplot(2,2,3)
        hold on; grid on;
        plot(t,Y(:,14))
        ylabel('dz/dk3')
        title('z variation wrt k3');
        subplot(2,2,4)
        hold on; grid on;
        plot(t,Y(:,19))
        ylabel('du/dk3')
        title('u variation wrt k3');
        suptitle('Vehavior of x,y,z,u variables wrt k3');
    
        figure(10)
        hold on; grid on;
        subplot(2,2,1)
        hold on; grid on;
        plot(t,Y(:,5))
        ylabel('dx/dk4')
        title('x variation wrt k4');
        subplot(2,2,2)
        hold on; grid on;
        plot(t,Y(:,10))
        ylabel('dy/dk4')
        title('y variation wrt k4');
        subplot(2,2,3)
        hold on; grid on;
        plot(t,Y(:,15))
        ylabel('dz/dk4')
        title('z variation wrt k4');
        subplot(2,2,4)
        hold on; grid on;
        plot(t,Y(:,20))
        ylabel('du/dk4')
        title('u variation wrt k4');
        suptitle('Behavior of x,y,z,u variables wrt k4');
    
       
end
 
% Add titles to figures in for loop1 above
figure(6)
hold on; grid on;
subplot(2,2,1)
hold on; grid on;
plot(t,Y(:,1))
ylabel('x')
title('x variation wrt varying Ks');
subplot(2,2,2)
hold on; grid on;
plot(t,Y(:,6))
ylabel('y')
title('y variation wrt varying Ks');
subplot(2,2,3)
hold on; grid on;
plot(t,Y(:,11))
ylabel('z')
title('z variation wrt varying Ks');
subplot(2,2,4)
hold on; grid on;
plot(t,Y(:,16))
ylabel('u')
title('u variation wrt varying Ks');
suptitle('Behavior of x,y,z,u variables for varying Ks');

figure(7)
hold on; grid on;
subplot(2,2,1)
hold on; grid on;
plot(t,Y(:,2))
ylabel('dx/dk1')
title('dx/dk1 variation wrt varying Ks');
subplot(2,2,2)
hold on; grid on;
plot(t,Y(:,7))
ylabel('dy/dk1')
title('dy/dk1 variation wrt varying Ks');
subplot(2,2,3)
hold on; grid on;
plot(t,Y(:,12))
ylabel('dz/dk1')
title('dz/dk1 variation wrt varying Ks');
subplot(2,2,4)
hold on; grid on;
plot(t,Y(:,17))
ylabel('du/dk1')
title('du/dk1 variation wrt varying Ks');
suptitle('Behavior of dx/dk1,dy/dk1,dz/dk1,du/dk1 for varying Ks');

figure(8)
hold on; grid on;
subplot(2,2,1)
hold on; grid on;
plot(t,Y(:,3))
ylabel('dx/dk2')
title('dx/dk2 variation wrt varying Ks');
subplot(2,2,2)
hold on; grid on;
plot(t,Y(:,8))
ylabel('dy/dk2')
title('dy/dk2 variation wrt varying Ks');
subplot(2,2,3)
hold on; grid on;
plot(t,Y(:,13))
ylabel('dz/dk2')
title('dz/dk2 variation wrt varying Ks');
subplot(2,2,4)
hold on; grid on;
plot(t,Y(:,18))
ylabel('du/dk2')
title('du/dk2 variation wrt varying Ks');
suptitle('Behavior of dx/dk2,dy/dk2,dz/dk2,du/dk2 for varying Ks');

figure(9)
hold on; grid on;
subplot(2,2,1)
hold on; grid on;
plot(t,Y(:,4))
ylabel('dx/dk3')
title('dx/dk3 variation wrt varying Ks');
subplot(2,2,2)
hold on; grid on;
plot(t,Y(:,9))
ylabel('dy/dk3')
title('dy/dk3 variation wrt varying Ks');
subplot(2,2,3)
hold on; grid on;
plot(t,Y(:,14))
ylabel('dz/dk3')
title('duz/dk3 variation wrt varying Ks');
subplot(2,2,4)
hold on; grid on;
plot(t,Y(:,19))
ylabel('du/dk3')
title('du/dk3 variation wrt varying Ks');
suptitle('Behavior of dx/dk3,dy/dk3,dz/dk3,du/dk3 for varying Ks');

figure(10)
hold on; grid on;
subplot(2,2,1)
hold on; grid on;
plot(t,Y(:,5))
ylabel('dx/dk4')
title('dx/dk4 variation wrt varying Ks');
subplot(2,2,2)
hold on; grid on;
plot(t,Y(:,10))
ylabel('dy/dk4')
title('dy/dk4 variation wrt varying Ks');
subplot(2,2,3)
hold on; grid on;
plot(t,Y(:,15))
ylabel('dz/dk4')
title('dz/dk4 variation wrt varying Ks');
subplot(2,2,4)
hold on; grid on;
plot(t,Y(:,20))
ylabel('du/dk4')
title('du/dk4 variation wrt varying Ks');
suptitle('Behavior of dx/dk4,dy/dk4,dz/dk4,du/dk4 for varying Ks');

h1 = figure(6);
h2 = figure(7);
h3 = figure(8);
h4 = figure(9);
h5 = figure(10);

saveas(h1,'Fig6_solnCurves_varyKs','jpg')
saveas(h2,'Fig7_variations2k1_varyKs','jpg')
saveas(h3,'Fig8_variations2k2_varyKs','jpg')
saveas(h4,'Fig9_variations2k3_varyKs','jpg')
saveas(h5,'Fig10_variations2k4_varyKs','jpg')

saveas(h1,'Fig6_solnCurves_varyKs','png')
saveas(h2,'Fig7_variations2k1_varyKs','png')
saveas(h3,'Fig8_variations2k2_varyKs','png')
saveas(h4,'Fig9_variations2k3_varyKs','png')
saveas(h5,'Fig10_variations2k4_varyKs','png')