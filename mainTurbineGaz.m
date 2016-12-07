function [state,Energy_losses,labels_Energy,etaMec,etaCyclen,etaToten,Exergy_losses,labels_Ex,etaRotex,etaCyclex,etaCombex,etaTotex,ma,mc,mg,lambda] = mainTurbineGaz(T1,r,etaPiC,kcc,T3,etaPiT,Pe,x,y,z,fuel)

T1 = T1 +273.15; % conversion °C --> K
T3 = T3 +273.15; % conversion °C --> K
Pe = Pe*10^3; % conversion MW --> kW

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% initialization %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

n = 4; % number of states
state = State_creation_Turbine(n);
% state table where 1 is the atmospheric air, 2 the compressor exit, 3 the
% exit of combustion chamber and 4 the exhaust gasses (after turbine)

%%%%%%%%%%%%%%
%%%% Data %%%%
%%%%%%%%%%%%%%

R = 8.3144621; %[J/(mol*K)]
kmec = 0.015;

%%%% pressure %%%%
state{1}.p = 1.01325; %[bar]
state{4}.p = 1.01325; %[bar]

%%%% air characteristics %%%%
% 21% O2 et 79% N2
Ma = 0.21*32 + 0.79*28; % [g/mol]
Ra = R/(Ma*10^(-3)); %[J/(kg*K)]

mO2_a = 0.21*(32/Ma); %[kg_O2/kg_air]
mN2_a = 0.79*(28/Ma); %[kg_N2/kg_air]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% State's calculation %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h0a = janaf('h','N2',273.15)*mN2_a + janaf('h','O2',273.15)*mO2_a;
s0a = janaf('s','N2',273.15)*mN2_a + janaf('s','O2',273.15)*mO2_a;

%%%% State 1 %%%%

state{1}.T = T1;
state{1}.h = janaf('h','N2',T1)*mN2_a + janaf('h','O2',T1)*mO2_a - h0a;
state{1}.s = janaf('s','N2',T1)*mN2_a + janaf('s','O2',T1)*mO2_a - s0a;

%%%% State 2 %%%%

state{2}.p = r*state{1}.p;

T2_0 = state{1}.T*(r)^((1.44-1)/1.44); % T1 en [K]
[T2,na] = solverCompressor(T2_0,etaPiC,Ra,T1,r,mO2_a,mN2_a);

state{2}.T = T2;
state{2}.h = janaf('h','N2',T2)*mN2_a + janaf('h','O2',T2)*mO2_a - h0a;
state{2}.s = janaf('s','N2',T2)*mN2_a + janaf('s','O2',T2)*mO2_a - s0a;

%%%% State 3 %%%%

state{3}.p = kcc*state{2}.p;
state{3}.T = T3;

%%%% States 3 & 4 - excess of air - mass flow %%%%
lambda0 = 2.3;
T4_0 = T3*(1/(r*kcc))^((1.44-1)/1.44); %[K]

% first iteration
[ng,T4] = solverTurbine(T4_0,etaPiT,state{3}.T,r,kcc,lambda0,x,y,z);

[mCO2_g,mH2O_g,mO2_g,mN2_g] = GasMassFraction(lambda0,x,y,z);

h0g = mCO2_g*janaf('h','CO2',273.15) + mH2O_g*janaf('h','H2O',273.15)...
    + mO2_g*janaf('h','O2',273.15) + mN2_g*janaf('h','N2',273.15);

state{3}.h = mCO2_g*janaf('h','CO2',T3) + mH2O_g*janaf('h','H2O',T3)...
    + mO2_g*janaf('h','O2',T3) + mN2_g*janaf('h','N2',T3) - h0g;

state{4}.h = mCO2_g*janaf('h','CO2',T4) + mH2O_g*janaf('h','H2O',T4)...
    + mO2_g*janaf('h','O2',T4) + mN2_g*janaf('h','N2',T4) - h0g;

[ma,mc,mg,lambda] = solverFlow(state,x,y,z,kmec,Pe,fuel);

lambda_prev = lambda0;
T4_prev = T4_0;

% iteration loop

while abs(lambda_prev - lambda)>0.01 || abs(T4 - T4_prev)>0.01
    lambda_prev = lambda;
    T4_prev = T4;
   
    [ng,T4] = solverTurbine(T4_0,etaPiT,state{3}.T,r,kcc,lambda,x,y,z);
    
    [mCO2_g,mH2O_g,mO2_g,mN2_g] = GasMassFraction(lambda,x,y,z);
    
    h0g = mCO2_g*janaf('h','CO2',273.15) + mH2O_g*janaf('h','H2O',273.15)...
    + mO2_g*janaf('h','O2',273.15) + mN2_g*janaf('h','N2',273.15);
    
    state{3}.h = mCO2_g*janaf('h','CO2',T3) + mH2O_g*janaf('h','H2O',T3)...
        + mO2_g*janaf('h','O2',T3) + mN2_g*janaf('h','N2',T3) - h0g;

    state{4}.h = mCO2_g*janaf('h','CO2',T4) + mH2O_g*janaf('h','H2O',T4)...
        + mO2_g*janaf('h','O2',T4) + mN2_g*janaf('h','N2',T4) - h0g;
   
    [ma,mc,mg,lambda] = solverFlow(state,x,y,z,kmec,Pe,fuel);
    
end

% States 3 and 4 
[mCO2_g,mH2O_g,mO2_g,mN2_g] = GasMassFraction(lambda,x,y,z);

h0g = mCO2_g*janaf('h','CO2',273.15) + mH2O_g*janaf('h','H2O',273.15)...
    + mO2_g*janaf('h','O2',273.15) + mN2_g*janaf('h','N2',273.15);

state{3}.h = mCO2_g*janaf('h','CO2',T3) + mH2O_g*janaf('h','H2O',T3)...
    + mO2_g*janaf('h','O2',T3) + mN2_g*janaf('h','N2',T3) - h0g;

state{4}.h = mCO2_g*janaf('h','CO2',T4) + mH2O_g*janaf('h','H2O',T4)...
    + mO2_g*janaf('h','O2',T4) + mN2_g*janaf('h','N2',T4) - h0g;


s0g = mCO2_g*janaf('s','CO2',273.15) + mH2O_g*janaf('s','H2O',273.15)...
    + mO2_g*janaf('s','O2',273.15) + mN2_g*janaf('s','N2',273.15);

state{3}.s = mCO2_g*janaf('s','CO2',T3) + mH2O_g*janaf('s','H2O',T3)...
    + mO2_g*janaf('s','O2',T3) + mN2_g*janaf('s','N2',T3) - s0g;

state{4}.s = mCO2_g*janaf('s','CO2',T4) + mH2O_g*janaf('s','H2O',T4)...
    + mO2_g*janaf('s','O2',T4) + mN2_g*janaf('s','N2',T4) - s0g;

state{4}.T = T4;

%%%% Exergy calculation %%%%
[state] = GasTurbineExergy(state);


%%%%%%%%%%%%%%%%%%
%%%% Analysis %%%%
%%%%%%%%%%%%%%%%%%

%%%% Energy analysis %%%%
[Energy_losses,labels_Energy,etaMec,etaCyclen,etaToten] = energyAnalysis4(state,ma,mc,mg,Pe,x,y,z,fuel);

%%%% Exergy analysis %%%%
[Exergy_losses,labels_Ex,etaMec,etaRotex,etaCyclex,etaCombex,etaTotex] = exergyAnalysis4(state,ma,mc,mg,x,y,z,Pe,fuel);





end

