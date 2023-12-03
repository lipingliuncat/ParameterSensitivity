function dxdt = SickleCell4DE(t,X)
% This Matlab code defines 4 DEs defining the dynamic
% of Sickle Cell model. 
% x'(t)=-k_1(x/C_1-1)z-k_2*x
% y'(t)=k_3(1-y/C_2)u+k_2*x
% z'(t)=k_1(x/C_1-1)z-k_4*z
% u'(t)=-k_3(1-y/C_2)u+k_4*z

global k1 k2 k3 k4 C_1 C_2
dxdt(1)= -k1*(X(1)/C_1-1)*X(3)-k2*X(1); 
dxdt(2)= k3*(1-(X(2)/C_2))*X(4)+k2*X(1);
dxdt(3)= k1*(X(1)/C_1-1)*X(3)-k4*X(3);
dxdt(4)= -k3*(1-(X(2)/C_2))*X(4)+k4*X(3);
dxdt = [dxdt(1),dxdt(2),dxdt(3),dxdt(4)]';
return