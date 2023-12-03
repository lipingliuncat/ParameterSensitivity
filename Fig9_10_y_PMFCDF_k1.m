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
    blockK = K_mat((b-1)*BlockSize+1:b * BlockSize,:);
    Allfobj((b-1)*BlockSize+1:b * BlockSize) = ...
        pMYRK4COSensMPSA(initX,blockK,t0,tf,dt,C_1,C_2,column,product);
end
%Store the sum of the total objective function
sum_fobj = sum(Allfobj);
% amin = min(Allfobj)
% amax = max(Allfobj)
% asize = size(Allfobj)
%Calculates the criterion value
%% Compute the criterion
criterion = prctile(Allfobj,[25,37.5,50,62.5,75]); %sum_fobj/n
criterion = [criterion,sum_fobj/n];
%% Plot out the CDF and PMF for the parameter
i=1;  % 1 for k_1, 2 for k_2, 3 for k_3 and 4 for k_4
for j=1:6
   criterion1 = criterion(j);
   mar_1 = K_mat(Allfobj <= criterion1,i);
   [f1,x1] =ecdf(mar_1);
   [n1,y1] = hist(mar_1);
   mar_2 = K_mat(Allfobj > criterion1,i);
   [f2,x2] =ecdf(mar_2);
   [n2,y2] = hist(mar_2);
   figure(1); 
   subplot(3,2,j); hold on; plot(x1,f1,'r-',x2,f2,'b-','LineWidth',1);
   figure(2); 
   subplot(3,2,j); hold on; plot(y1,n1/1000,'rx-',y2,n2/1000,'bx-','LineWidth',1);
end

figure(1); 
subplot(3,2,1); box on; grid on; axis([0.01 0.06 0 1]);
xlabel('$k_1$','Interpreter','latex');
ylabel('CDF'); txt={'(a)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);
subplot(3,2,2); box on; grid on; axis([0.01 0.06 0 1]);
xlabel('$k_1$','Interpreter','latex');
ylabel('CDF'); txt={'(b)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);
subplot(3,2,3); box on; grid on; axis([0.01 0.06 0 1]);
xlabel('$k_1$','Interpreter','latex');
ylabel('CDF'); txt={'(c)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);
subplot(3,2,4); box on; grid on; axis([0.01 0.06 0 1]);
xlabel('$k_1$','Interpreter','latex');
ylabel('CDF'); txt={'(d)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);
subplot(3,2,5); box on; grid on; axis([0.01 0.06 0 1]);
xlabel('$k_1$','Interpreter','latex');
ylabel('CDF'); txt={'(e)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);
subplot(3,2,6); box on; grid on; axis([0.01 0.06 0 1]);
xlabel('$k_1$','Interpreter','latex');
ylabel('CDF'); txt={'(f)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);

figure(2); 
subplot(3,2,1); box on; grid on; axis([0.01 0.06 0 1]);
xlabel('$k_1$','Interpreter','latex');
ylabel('PMF'); txt={'(a)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);
subplot(3,2,2); box on; grid on; axis([0.01 0.06 0 1]);
xlabel('$k_1$','Interpreter','latex');
ylabel('PMF'); txt={'(b)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);
subplot(3,2,3); box on; grid on; axis([0.01 0.06 0 1]);
xlabel('$k_1$','Interpreter','latex');
ylabel('PMF'); txt={'(c)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);
subplot(3,2,4); box on; grid on; axis([0.01 0.06 0 1]);
xlabel('$k_1$','Interpreter','latex');
ylabel('PMF'); txt={'(d)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);
subplot(3,2,5); box on; grid on; axis([0.01 0.06 0 1]);
xlabel('$k_1$','Interpreter','latex');
ylabel('PMF'); txt={'(e)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);
subplot(3,2,6); box on; grid on; axis([0.01 0.06 0 1]);
xlabel('$k_1$','Interpreter','latex');
ylabel('PMF'); txt={'(f)'}; text(max(xlim)*0.2,max(ylim)*0.9,txt);

ax1 = figure(1);
ax2 = figure(2);
exportgraphics(ax1,'Fig10.png','Resolution',1200);
exportgraphics(ax2,'Fig9.png','Resolution',1200);