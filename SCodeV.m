function dxdt = SCodeV(X,K,C_1,C_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code defines an extended set of 20 PDEs defining the dynamic
% of Sickle Cell model and the variations of the dependant variables with
% respect to four parameters k1,k2,k3 and k4. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = size(K,1);
    dxdt = zeros(n,4);
    T1 =  K(:,2).*X(:,1);
    T2 =  K(:,4).*X(:,3);
    T3 =  K(:,1).*(X(:,1)/C_1-1).*X(:,3);
    T4 = -K(:,3).*(X(:,2)/C_2-1).*X(:,4);
    
    dxdt(:,1)= -T3 - T1; 
    dxdt(:,2)=  T4 + T1;
    dxdt(:,3)=  T3 - T2;
    dxdt(:,4)= -T4 + T2;
end