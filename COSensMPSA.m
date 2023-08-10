function dxdt = COSensMPSA(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code defines an extended set of 20 PDEs defining the dynamic
% of Sickle Cell model and the variations of the dependant variables with
% respect to four parameters k1,k2,k3 and k4. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global k1 k2 k3 k4 C_1 C_2

dxdt(1)= -k1*(X(1)/C_1-1)*X(3)-k2*X(1); 

dxdt(2)= k3*(1-(X(2)/C_2))*X(4)+k2*X(1);


dxdt(3)= k1*(X(1)/C_1-1)*X(3)-k4*X(3);


dxdt(16)= -k3*(1-(X(2)/C_2))*X(4)+k4*X(3);


dxdt = [dxdt(1),dxdt(2),dxdt(3),dxdt(4)];
end