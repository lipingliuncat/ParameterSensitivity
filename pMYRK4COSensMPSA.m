function difference = pMYRK4COSensMPSA(X,K,t0,tf,dt,C_1,C_2,column,product)
n_max = floor((tf-t0)/dt);
difference = (X(:,column)-product(1)).^2;%/1000;
t = zeros(n_max+1,1); t(1) = t0;
for n = 1:1:n_max    
    K1 = pCOSensMPSA(X,K,C_1,C_2);
    K2 = pCOSensMPSA(X+dt/2*K1,K,C_1,C_2);
    K3 = pCOSensMPSA(X+dt/2*K2,K,C_1,C_2);
    K4 = pCOSensMPSA(X+dt*K3,K,C_1,C_2);  
    X = X + (K1+2*K2+2*K3+K4)*(dt/6);
    difference = difference + ((X(:,column)-product(n+1)).^2);%/1000;
    t(n+1) = n * dt;
end
end
