function [ng,T4] = solverTurbine(T4_0,etaPiT,T3,r,kcc,lambda,x,y,z)


R = 8.3144621; %[J/(mol*K)]

% Gas'components mass fraction
[mCO2_g,mH2O_g,mO2_g,mN2_g,Mg] = GasMassFraction(lambda,x,y,z);

Rg = R/(Mg*10^(-3));

%%%% Iterative solver %%%%

% step 1
T = linspace(T4_0,T3,500);
cp = mean(mCO2_g*janaf('c','CO2',T) + mH2O_g*janaf('c','H2O',T)...
    + mO2_g*janaf('c','O2',T) + mN2_g*janaf('c','N2',T))*10^3;

ng = (1-(Rg*etaPiT/cp))^(-1);
T4 = T3*(1/(kcc*r))^((ng-1)/ng);

T4_prev = T4_0;
ng_prev = ng +5;

% while loop
while (abs(T4-T4_prev) > 0.0001) || (abs(ng-ng_prev) > 0.0001)
    T4_prev = T4;
    ng_prev = ng;
    
    T = linspace(T4,T3,500);
    cp = mean(mCO2_g*janaf('c','CO2',T) + mH2O_g*janaf('c','H2O',T)...
        + mO2_g*janaf('c','O2',T) + mN2_g*janaf('c','N2',T))*10^3;

    ng = (1-(Rg*etaPiT/cp))^(-1);
    T4 = T3*(1/(kcc*r))^((ng-1)/ng);
end