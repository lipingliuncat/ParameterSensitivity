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
subplot(2,2,1); hold on; grid on; plot(t,Y(:,1)-Y0(:,1));
subplot(2,2,2); hold on; grid on; plot(t,Y(:,6)-Y0(:,6));
subplot(2,2,3); hold on; grid on; plot(t,Y(:,11)-Y0(:,11));
subplot(2,2,4); hold on; grid on; plot(t,Y(:,16)-Y0(:,16));
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
    [t,Y] = ode45(@SickleCellDE,tspan,X0);
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
    figure(1); 
    subplot(2,2,1);     plot(t,Y(:,1)-Y0(:,1));
    subplot(2,2,2);     plot(t,Y(:,6)-Y0(:,6));
    subplot(2,2,3);     plot(t,Y(:,11)-Y0(:,11));
    subplot(2,2,4);     plot(t,Y(:,16)-Y0(:,16));
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
figure(1); 
subplot(2,2,1); ylabel('x-x_{ref}');
title('x variation wrt varying Ks');
subplot(2,2,2); ylabel('y-y_{ref}'); 
title('y variation wrt varying Ks');
subplot(2,2,3); ylabel('z-z_{ref}');
title('z variation wrt varying Ks');
subplot(2,2,4); ylabel('u-u_{ref}');
title('u variation wrt varying Ks');
sgtitle('Behavior of x,y,z,u variables for varying Ks');
h1 = figure(1);
RefSolnDiff_range = [x_diff_min,y_diff_min,z_diff_min,u_diff_min;...
    x_diff_max,y_diff_max,z_diff_max,u_diff_max;...
    x_diff_max-x_diff_min,y_diff_max-y_diff_min,z_diff_max-z_diff_min,...
    u_diff_max-u_diff_min];
fileID = fopen('SolnDiffRange.txt','w');
fprintf(fileID,'%10s %12s %12s %12s\r\n','x_diff','y_diff',...
    'z_diff','u_diff');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\r\n',RefSolnDiff_range');
fclose(fileID);
k1_range = [k1x_min,k1y_min,k1z_min,k1u_min;
    k1x_max,k1y_max,k1z_max,k1u_max;
    k1x_max-k1x_min,k1y_max-k1y_min,k1z_max-k1z_min,k1u_max-k1u_min];
fileID = fopen('K1Range.txt','w');
fprintf(fileID,'%10s %12s %12s %12s\r\n','k1x','k1y','k1z','k1u');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\r\n',k1_range');
fclose(fileID);
k2_range = [k2x_min,k2y_min,k2z_min,k2u_min;
    k2x_max,k2y_max,k2z_max,k2u_max;
    k2x_max-k2x_min,k2y_max-k2y_min,k2z_max-k2z_min,k2u_max-k2u_min];
fileID = fopen('K2Range.txt','w');
fprintf(fileID,'%10s %12s %12s %12s\r\n','k2x','k2y','k2z','k2u');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\r\n',k2_range');
fclose(fileID);
k3_range = [k3x_min,k3y_min,k3z_min,k3u_min;
    k3x_max,k3y_max,k3z_max,k3u_max;
    k3x_max-k3x_min,k3y_max-k3y_min,k3z_max-k3z_min,k3u_max-k3u_min];
fileID = fopen('K3Range.txt','w');
fprintf(fileID,'%10s %12s %12s %12s\r\n','k3x','k3y','k3z','k3u');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\r\n',k3_range');
fclose(fileID);
k4_range = [k4x_min,k4y_min,k4z_min,k4u_min;
    k4x_max,k4y_max,k4z_max,k4u_max;
    k4x_max-k4x_min,k4y_max-k4y_min,k4z_max-k4z_min,k4u_max-k4u_min];
fileID = fopen('K4Range.txt','w');
fprintf(fileID,'%10s %12s %12s %12s\r\n','k4x','k4y','k4z','k4u');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\r\n',k4_range');
fclose(fileID);
% For loop1: sensitivity of x,y,z,u wrt k1,k2,k3,k4 when all four 
% parametersvary simultaneously
figure(7); 
subplot(2,2,1); hold on; grid on;     
ylabel('x');
title('x variation wrt varying Ks');
subplot(2,2,2); hold on; grid on;  
ylabel('y');
title('y variation wrt varying Ks');
subplot(2,2,3); hold on; grid on;   
ylabel('z');
title('z variation wrt varying Ks');
subplot(2,2,4); hold on; grid on;   
ylabel('u');
title('u variation wrt varying Ks');
sgtitle('Behavior of x,y,z,u variables for varying Ks');
figure(8); 
subplot(2,2,1); hold on; grid on; 
ylabel('dx/dk_1')
title('dx/dk_1 variation wrt varying Ks');
subplot(2,2,2); hold on; grid on;
ylabel('dy/dk_1');
title('dy/dk_1 variation wrt varying Ks');
subplot(2,2,3); hold on; grid on;
ylabel('dz/dk_1');
title('dz/dk_1 variation wrt varying Ks');
subplot(2,2,4); hold on; grid on;
ylabel('du/dk_1');
title('du/dk_1 variation wrt varying Ks');
sgtitle('Behavior of dx/dk_1,dy/dk_1,dz/dk_1,du/dk_1 for varying Ks');
figure(9); 
subplot(2,2,1); hold on; grid on;
ylabel('dx/dk_2');
title('dx/dk_2 variation wrt varying Ks');
subplot(2,2,2); hold on; grid on;
ylabel('dy/dk_2');
title('dy/dk_2 variation wrt varying Ks');
subplot(2,2,3); hold on; grid on;
ylabel('dz/dk_2');
title('dz/dk_2 variation wrt varying Ks');
subplot(2,2,4); hold on; grid on;
ylabel('du/dk_2');
title('du/dk_2 variation wrt varying Ks');
sgtitle('Behavior of dx/dk_2,dy/dk_2,dz/dk_2,du/dk_2 for varying Ks');
figure(10); 
subplot(2,2,1); hold on; grid on;
ylabel('dx/dk_3'); 
title('dx/dk_3 variation wrt varying Ks');
subplot(2,2,2); hold on; grid on;
ylabel('dy/dk_3');
title('dy/dk_3 variation wrt varying Ks');
subplot(2,2,3); hold on; grid on;
ylabel('dz/dk_3');
title('dz/dk_3 variation wrt varying Ks');
subplot(2,2,4); hold on; grid on;
ylabel('du/dk_3'); 
title('du/dk_3 variation wrt varying Ks');
sgtitle('Behavior of dx/dk_3,dy/dk_3,dz/dk_3,du/dk_3 for varying Ks');
figure(11); 
subplot(2,2,1); hold on; grid on;
ylabel('dx/dk_4');
title('dx/dk_4 variation wrt varying Ks');
subplot(2,2,2); hold on; grid on;
ylabel('dy/dk_4');
title('dy/dk_4 variation wrt varying Ks');
subplot(2,2,3); hold on; grid on;
ylabel('dz/dk_4');
title('dz/dk_4 variation wrt varying Ks');
subplot(2,2,4); hold on; grid on;
ylabel('du/dk_4'); 
title('du/dk_4 variation wrt varying Ks');
sgtitle('Behavior of dx/dk_4,dy/dk_4,dz/dk_4,du/dk_4 for varying Ks');

for i_k = 1:1:dim
        k1 = k1v(i_k);
        k2 = k2v(i_k);
        k3 = k3v(i_k);
        k4 = k4v(i_k);
        [t,Y] = ode45(@SickleCellDE,tspan,X0);
        figure(7); 
        subplot(2,2,1); plot(t,Y(:,1));
        subplot(2,2,2); plot(t,Y(:,6));
        subplot(2,2,3); plot(t,Y(:,11)); 
        subplot(2,2,4); plot(t,Y(:,16)); 
        figure(8); 
        subplot(2,2,1); plot(t,Y(:,2)); 
        subplot(2,2,2); plot(t,Y(:,7)); 
        subplot(2,2,3); plot(t,Y(:,12)); 
        subplot(2,2,4); plot(t,Y(:,17)); 
        figure(9);      
        subplot(2,2,1); plot(t,Y(:,3)); 
        subplot(2,2,2); plot(t,Y(:,8));
        subplot(2,2,3); plot(t,Y(:,13)); 
        subplot(2,2,4); plot(t,Y(:,18));
        figure(10);     
        subplot(2,2,1); plot(t,Y(:,4));
        subplot(2,2,2); plot(t,Y(:,9));
        subplot(2,2,3); plot(t,Y(:,14));
        subplot(2,2,4); plot(t,Y(:,19));
        figure(11);     
        subplot(2,2,1); plot(t,Y(:,5));
        subplot(2,2,2); plot(t,Y(:,10));
        subplot(2,2,3); plot(t,Y(:,15));
        subplot(2,2,4); plot(t,Y(:,20));
end

h1 = figure(7);
h2 = figure(8);
h3 = figure(9);
h4 = figure(10);
h5 = figure(11);
saveas(h1,'Fig7_solnCurves_varyKs','jpg');
saveas(h2,'Fig8_variations2k1_varyKs','jpg');
saveas(h3,'Fig9_variations2k2_varyKs','jpg');
saveas(h4,'Fig10_variations2k3_varyKs','jpg');
saveas(h5,'Fig11_variations2k4_varyKs','jpg');
saveas(h1,'Fig7_solnCurves_varyKs','png');
saveas(h2,'Fig8_variations2k1_varyKs','png');
saveas(h3,'Fig9_variations2k2_varyKs','png');
saveas(h4,'Fig10_variations2k3_varyKs','png');
saveas(h5,'Fig11_variations2k4_varyKs','png');