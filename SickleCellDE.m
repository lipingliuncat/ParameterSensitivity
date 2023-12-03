function dxdt = SickleCellDEr(t,X)
% This Matlab code defines 4 DEs defining the dynamic
% of Sickle Cell model. 
global k1 k2 k3 k4 C_1 C_2
xt=X(1); xk1=X(2); xk2=X(3); xk3=X(4); xk4=X(5);
yt=X(6); yk1=X(7); yk2=X(8); yk3=X(9); yk4=X(10);
zt=X(11); zk1=X(12); zk2=X(13); zk3=X(14); zk4=X(15);
ut=X(16); uk1=X(17); uk2=X(18); uk3=X(19); uk4=X(20);
dxdt(1)= -k1*(xt/C_1-1)*zt-k2*xt;
dxdt(2)=(1-xt/C_1)*zt+k1*(1-xt/C_1)*zk1-k1/C_1*xk1*zt-k2*xk1;
dxdt(3)=k1*(1-xt/C_1)*zk2-k1/C_1*xk2*zt-k2*xk2-xt;
dxdt(4)=k1*(1-xt/C_1)*zk3-k1/C_1*xk3*zt-k2*xk3;
dxdt(5)=k1*(1-xt/C_1)*zk4-k1/C_1*xk4*zt-k2*xk4;
dxdt(6)= k3*(1-yt/C_2)*ut+k2*xt;
dxdt(7)=k3*(1-yt/C_2)*uk1-k3/C_2*yk1*ut+k2*xk1;
dxdt(8)=k3*(1-yt/C_2)*uk2-k3/C_2*yk2*ut+k2*xk2+xt;
dxdt(9)=(1-yt/C_2)*ut+k3*(1-yt/C_2)*uk3-k3/C_2*yk3*ut+k2*xk3;
dxdt(10)=k3*(1-yt/C_2)*uk4-k3/C_2*yk4*ut+k2*xk4;
dxdt(11)= k1*(xt/C_1-1)*zt-k4*zt;
dxdt(12)=(xt/C_1-1)*zt+k1*(xt/C_1-1)*zk1+k1/C_1*xk1*zt-k4*zk1;
dxdt(13)=k1*(xt/C_1-1)*zk2+k1/C_1*xk2*zt-k4*zk2;
dxdt(14)=k1*(xt/C_1-1)*zk3+k1/C_1*xk3*zt-k4*zk3;
dxdt(15)=k1*(xt/C_1-1)*zk4+k1/C_1*xk4*zt-k4*zk4-zt;
dxdt(16)= -k3*(1-yt/C_2)*ut+k4*zt;
dxdt(17)= k3/C_2*yk1*ut+k3*(yt/C_2-1)*uk1+k4*zk1;
dxdt(18)= -k3*(1-yt/C_2)*uk2+k3/C_2*yk2*ut+k4*zk2;
dxdt(19)= -(1-yt/C_2)*ut-k3*(1-yt/C_2)*uk3+k3/C_2*yk3*ut+k4*zk3;
dxdt(20)= -k3*(1-yt/C_2)*uk4+k3/C_2*yk4*ut+k4*zk4+zt;
dxdt = [dxdt(1),dxdt(2),dxdt(3),dxdt(4),dxdt(5),dxdt(6),dxdt(7),dxdt(8),...
    dxdt(9),dxdt(10),dxdt(11),dxdt(12),dxdt(13),dxdt(14),dxdt(15),...
    dxdt(16),dxdt(17),dxdt(18),dxdt(19),dxdt(20)]';
return