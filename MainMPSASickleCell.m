%Initialization of the parameter of the model
k1=.028;
k2=.07;
k3=.1;
k4=.01;
K = [k1 k2 k3 k4];
C_1=.4;
C_2=.1;
t0=0;
tf=250;
dt=0.01;
%Calling My runge-kutta 4 method to compute the original (nominal)
%behaviour
X0=[0.0036,0.0,.1750,0];
[t,Y] =MYRK4COSensMPSA(X0,K,t0,tf,dt,C_1,C_2);
%To compute “x” column must = 1, “y” column must = 2, “z” column 
%must = 3 and “u” column must equal 4.
column = 2;  % The considered output is x, the 1st variable
product=squeeze(Y(:,column,:))';
figure(9);
hold on; grid on;
plot(t,squeeze(Y(:,1,:)),'b')
plot(t,squeeze(Y(:,2,:)),'r')
plot(t,squeeze(Y(:,3,:)),'g')
plot(t,squeeze(Y(:,4,:)),'k')
legend('De-oxy Monomers[x(1)]','CO-bound Monomers [x(2)]',...
'De-oxy polymers [x(3)]','CO-bound Polymers [x(4)]');
title('Molar Cencentration') % put a title on the plot
xlabel('Time [min]')%  label the x-axis
ylabel('Population [Mm]') %  label the y-axis   
%is the number of the point we would like to use for the
%simulation
%% n = 1000000;  original setting
n = 10000;
K_init = 0.5 * K;
K_fin = 2 * K;
K_mat = zeros(n,length(K));
%Creates a uniform random number distribution between the set of 2 points
for i = 1:length(K)
    K_mat(:,i) = random('unif',K_init(i),K_fin(i),n,1);
end
%Stores the objective function of the simulation
Allfobj=zeros(n,1);
%Simulation
%% BlockSize = 10000; original setting
BlockSize = 1000;
Blocks = n / BlockSize;
initX = repmat([0.0036,0.0,0.1750,0.0],BlockSize,1);
for b = 1:Blocks
    b
    tic
    blockK = K_mat((b-1)*BlockSize+1:b * BlockSize,:);
    Allfobj((b-1)*BlockSize+1:b * BlockSize) = ...
        pMYRK4COSensMPSA(initX,blockK,t0,tf,dt,C_1,C_2,column,product);
    toc
end
%Store the sum of the total objective function
sum_fobj = sum(Allfobj);
% amin = min(Allfobj)
% amax = max(Allfobj)
% asize = size(Allfobj)
save data_5;
%Calculates the criterion value
%% Compute the criterion
criterion = prctile(Allfobj,[25,37.5,50,62.5,75]); %sum_fobj/n
criterion = [criterion,sum_fobj/n];
%% Plot out the CDF for four parameters
h = [];
p = [];
rank = [];
for i=1:1:4
    h2 = [];
    p2 = [];
    covk2 = [];
    rank2 = [];
    for j=1:6
        criterion1 = criterion(j);
        mar_1 = K_mat(Allfobj <= criterion1,i);
        [f1,x1] =ecdf(mar_1);
        [n1,y1] = hist(mar_1);
        mar_2 = K_mat(Allfobj > criterion1,i);
        [f2,x2] =ecdf(mar_2);
        [n2,y2] = hist(mar_2);
        [h1,p1,rank1]= kstest2(mar_1,mar_2,0.05);
        h2 = [h2;h1];
        p2 = [p2;p1];
        rank2 = [rank2;rank1];
        figure(i); 
        hold on; grid on;
        subplot(3,2,j)
        hold on; grid on;
        plot(x1,f1,'r-',x2,f2,'b-');
        %legend('Acceptable','Unacceptable');
        figure(i+4)
        hold on; grid on;
        subplot(3,2,j)
        hold on; grid on;
        plot(y1,n1/1000,'rx-',y2,n2/1000,'bx-');
        %legend('Acceptable (Red)','Unacceptable(Blue)');
    end
    h = [h,h2];
    p = [p,p2];
    rank = [rank,rank2];
end 
figure(1); 
subplot(3,2,1)
hold on; grid on;
xlabel('Parameter k1'); 
ylabel('25-prctile Obj ECDFs');
subplot(3,2,2)
hold on; grid on;
xlabel('Parameter k1'); 
ylabel('37.5-prctile Obj ECDFs');
subplot(3,2,3)
hold on; grid on;
xlabel('Parameter k1'); 
ylabel('50-prctile Obj ECDFs');
subplot(3,2,4)
hold on; grid on;
xlabel('Parameter k1'); 
ylabel('62.5-prctile Obj ECDFs');
subplot(3,2,5)
hold on; grid on;
xlabel('Parameter k1'); 
ylabel('75-prctile ObjECDFs');
subplot(3,2,6)
hold on; grid on;
xlabel('Parameter k1'); 
ylabel('Mean Obj ECDFs');
sgtitle('Acceptable (Red) and Unacceptable (Blue) ECDFs of k1');
figure(2); 
subplot(3,2,1)
hold on; grid on;
xlabel('Parameter k2'); 
ylabel('25-prctile Obj ECDFs');
subplot(3,2,2)
hold on; grid on;
xlabel('Parameter k2'); 
ylabel('37.5-prctile Obj ECDFs');
subplot(3,2,3)
hold on; grid on;
xlabel('Parameter k2'); 
ylabel('50-prctile Obj ECDFs');
subplot(3,2,4)
hold on; grid on;
xlabel('Parameter k2'); 
ylabel('62.5-prctile Obj ECDFs');
subplot(3,2,5)
hold on; grid on;
xlabel('Parameter k2'); 
ylabel('75-prctile ObjE CDFs');
subplot(3,2,6)
hold on; grid on;
xlabel('Parameter k2'); 
ylabel('Mean Obj ECDFs');
sgtitle('Acceptable (Red) and Unacceptable (Blue) ECDFs of k2');
figure(3); 
subplot(3,2,1)
hold on; grid on;
xlabel('Parameter k3'); 
ylabel('25-prctile Obj ECDFs');
subplot(3,2,2)
hold on; grid on;
xlabel('Parameter k3'); 
ylabel('37.5-prctile Obj ECDFs');
subplot(3,2,3)
hold on; grid on;
xlabel('Parameter k3'); 
ylabel('50-prctile Obj ECDFs');
subplot(3,2,4)
hold on; grid on;
xlabel('Parameter k3'); 
ylabel('62.5-prctile Obj ECDFs');
subplot(3,2,5)
hold on; grid on;
xlabel('Parameter k3'); 
ylabel('75-prctile Obj ECDFs');
subplot(3,2,6)
hold on; grid on;
xlabel('Parameter k3'); 
ylabel('Mean Obj ECDFs');
sgtitle('Acceptable (Red) and Unacceptable (Blue) ECDFs of k3');
figure(4); 
subplot(3,2,1)
hold on; grid on;
xlabel('Parameter k4'); 
ylabel('25-prctile Obj ECDFs');
subplot(3,2,2)
hold on; grid on;
xlabel('Parameter k4'); 
ylabel('37.5-prctile Obj ECDFs');
subplot(3,2,3)
hold on; grid on;
xlabel('Parameter k4'); 
ylabel('50-prctile Obj ECDFs');
subplot(3,2,4)
hold on; grid on;
xlabel('Parameter k4'); 
ylabel('62.5-prctile Obj ECDFs');
subplot(3,2,5)
hold on; grid on;
xlabel('Parameter k4'); 
ylabel('75-prctile Obj ECDFs');
subplot(3,2,6)
hold on; grid on;
xlabel('Parameter k4'); 
ylabel('Mean Obj ECDFs');
sgtitle('Acceptable (Red) and Unacceptable (Blue) ECDFs of k4');
figure(5);
subplot(3,2,1)
hold on; grid on;
xlabel('Parameter k1'); 
ylabel('25-prctile Obj PMFs');
subplot(3,2,2)
hold on; grid on;
xlabel('Parameter k1'); 
ylabel('37.5-prctile Obj PMFs');
subplot(3,2,3)
hold on; grid on;
xlabel('Parameter k1'); 
ylabel('50-prctile Obj PMFs');
subplot(3,2,4)
hold on; grid on;
xlabel('Parameter k1'); 
ylabel('62.5-prctile Obj PMFs');
subplot(3,2,5)
hold on; grid on;
xlabel('Parameter k1'); 
ylabel('75-prctile Obj PMFs');
subplot(3,2,6)
hold on; grid on;
xlabel('Parameter k1'); 
ylabel('Mean Obj PMFs');
sgtitle('Acceptable (Red) and Unacceptable (Blue) PMFs of k1');
figure(6); 
subplot(3,2,1)
hold on; grid on;
xlabel('Parameter k2'); 
ylabel('25-prctile Obj PMFs');
subplot(3,2,2)
hold on; grid on;
xlabel('Parameter k2'); 
ylabel('37.5-prctile Obj PMFs');
subplot(3,2,3)
hold on; grid on;
xlabel('Parameter k2'); 
ylabel('50-prctile Obj PMFs');
subplot(3,2,4)
hold on; grid on;
xlabel('Parameter k2'); 
ylabel('62.5-prctile Obj PMFs');
subplot(3,2,5)
hold on; grid on;
xlabel('Parameter k2'); 
ylabel('75-prctile Obj PMFs');
subplot(3,2,6)
hold on; grid on;
xlabel('Parameter k2'); 
ylabel('Mean Obj PMFs');
sgtitle('Acceptable (Red) and Unacceptable (Blue) PMFs of k2');
figure(7); 
subplot(3,2,1)
hold on; grid on;
xlabel('Parameter k3'); 
ylabel('25-prctile Obj PMFs');
subplot(3,2,2)
hold on; grid on;
xlabel('Parameter k3'); 
ylabel('37.5-prctile Obj PMFs');
subplot(3,2,3)
hold on; grid on;
xlabel('Parameter k3'); 
ylabel('50-prctile Obj PMFs');
subplot(3,2,4)
hold on; grid on;
xlabel('Parameter k3'); 
ylabel('62.5-prctile Obj PMFs');
subplot(3,2,5)
hold on; grid on;
xlabel('Parameter k3'); 
ylabel('75-prctile Obj PMFs');
subplot(3,2,6)
hold on; grid on;
xlabel('Parameter k3'); 
ylabel('Mean Obj PMFs');
sgtitle('Acceptable (Red) and Unacceptable (Blue) PMFs of k3');
figure(8); 
subplot(3,2,1)
hold on; grid on;
xlabel('Parameter k4'); 
ylabel('25-prctile Obj PMFs');
subplot(3,2,2)
hold on; grid on;
xlabel('Parameter k4'); 
ylabel('37.5-prctile Obj PMFs');
subplot(3,2,3)
hold on; grid on;
xlabel('Parameter k4'); 
ylabel('50-prctile Obj PMFs');
subplot(3,2,4)
hold on; grid on;
xlabel('Parameter k4'); 
ylabel('62.5-prctile Obj PMFs');
subplot(3,2,5)
hold on; grid on;
xlabel('Parameter k4'); 
ylabel('75-prctile Obj PMFs');
subplot(3,2,6)
hold on; grid on;
xlabel('Parameter k4'); 
ylabel('Mean Obj PMFs');
sgtitle('Acceptable (Red) and Unacceptable (Blue) PMFs of k4');
h1 = figure(1);
h2 = figure(2);
h3 = figure(3);
h4 = figure(4);
h5 = figure(5);
h6 = figure(6);
h7 = figure(7);
h8 = figure(8);
saveas(h1,'AccepRejectECDF4K1','jpg');
saveas(h2,'AccepRejectECDF4K2','jpg');
saveas(h3,'AccepRejectECDF4K3','jpg');
saveas(h4,'AccepRejectECDF4K4','jpg');
saveas(h5,'AccepRejectPMF4K1','jpg');
saveas(h6,'AccepRejectPMF4K2','jpg');
saveas(h7,'AccepRejectPMF4K3','jpg');
saveas(h8,'AccepRejectPMF4K4','jpg');
display('The sensitivity ranks at 5% significance level are:');
rank
fileID = fopen('SensivityRank.txt','w');
fprintf(fileID,'%10s %12s %12s %12s\r\n','k1','k2',...
    'k3','k4');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\r\n',rank');
fclose(fileID);
display('The K-S test outcomes at 5% significance level are:');
h
fileID = fopen('KSTestResult.txt','w');
fprintf(fileID,'%10s %12s %12s %12s\r\n','Null or Alt','Null or Alt',...
    'Null or Alt','Null or Alt');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\r\n',h');
fclose(fileID);
display('The K-S test asymptotic p-values at 5% significance level are:');
p
fileID = fopen('pValue.txt','w');
fprintf(fileID,'%10s %12s %12s %12s\r\n','p-vlaue','p-vlaue',...
    'p-vlaue','p-vlaue');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\r\n',p');
fclose(fileID);
