function [state,Energy_losses,labels_Energy,etaMec,etaCyclen,etaToten,Exergy_losses,labels_Ex,etaRotex,etaCyclex,etaCombex,etaTotex,ma,mc,mg,lambda] = mainTurbineGazRecup(T1,r,etaPiC,kcc,T4,etaPiT,Pe,x,y,z,NTU)

T1 = T1 +273.15; % conversion °C --> K
T4 = T4 +273.15; % conversion °C --> K
Pe = Pe*10^3; % conversion MW --> kW

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% initialization %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

n = 6; % number of states
state = State_creation_Turbine(n);
% state table where 1 is the atmospheric air, 2 the compressor exit, 3 the air after recuperator, 4 the
% exit of combustion chamber, 5 the gas after recuperator and 6 the exhaust gasses (after turbine)

%%%%%%%%%%%%%%
%%%% Data %%%%
%%%%%%%%%%%%%%

R = 8.3144621; %[J/(mol*K)]
kmec = 0.015;

%%%% pressure %%%%
state{1}.p = 1.01325; %[bar]
state{6}.p = 1.01325; %[bar]

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

%%%% State 3 (2R) %%%%

state{3}.p = 0.95*state{2}.p;

%%%% State 3,4,5 - mass flows - air excess %%%%
state{4}.T = T4;
state{4}.p = state{3}.p*kcc;
state{5}.p = state{6}.p/0.95;

lambda0 = 2.3;
T5_0 = T4*(1/(0.95*r*kcc))^((1.44-1)/1.44); %[K]

% first iteration
[ng,T5] = solverTurbineRecup(T5_0,etaPiT,state{4}.T,r,kcc,lambda0,x,y,z);

[mCO2_g,mH2O_g,mO2_g,mN2_g] = GasMassFraction(lambda0,x,y,z);

h0g = mCO2_g*janaf('h','CO2',273.15) + mH2O_g*janaf('h','H2O',273.15)...
    + mO2_g*janaf('h','O2',273.15) + mN2_g*janaf('h','N2',273.15);

state{4}.h = mCO2_g*janaf('h','CO2',T4) + mH2O_g*janaf('h','H2O',T4)...
    + mO2_g*janaf('h','O2',T4) + mN2_g*janaf('h','N2',T4) - h0g;

state{5}.h = mCO2_g*janaf('h','CO2',T5) + mH2O_g*janaf('h','H2O',T5)...
    + mO2_g*janaf('h','O2',T5) + mN2_g*janaf('h','N2',T5) - h0g;

[ma,mc,mg,lambda,T3,state] = solverFlowRecup(state,x,y,z,kmec,Pe,NTU,T5);

lambda_prev = lambda0;
T5_prev = T5_0;

% iteration loop


while abs(lambda_prev - lambda)>0.01 || abs(T5 - T5_prev)>0.01
    lambda_prev = lambda;
    T5_prev = T5;
   
    [ng,T5] = solverTurbineRecup(T5,etaPiT,state{4}.T,r,kcc,lambda,x,y,z);
    [mCO2_g,mH2O_g,mO2_g,mN2_g] = GasMassFraction(lambda,x,y,z);
    
    h0g = mCO2_g*janaf('h','CO2',273.15) + mH2O_g*janaf('h','H2O',273.15)...
    + mO2_g*janaf('h','O2',273.15) + mN2_g*janaf('h','N2',273.15);
    
    state{4}.h = mCO2_g*janaf('h','CO2',T4) + mH2O_g*janaf('h','H2O',T4)...
        + mO2_g*janaf('h','O2',T4) + mN2_g*janaf('h','N2',T4) - h0g;

    state{5}.h = mCO2_g*janaf('h','CO2',T5) + mH2O_g*janaf('h','H2O',T5)...
        + mO2_g*janaf('h','O2',T5) + mN2_g*janaf('h','N2',T5) - h0g;
    
    [ma,mc,mg,lambda,T3,state] = solverFlowRecup(state,x,y,z,kmec,Pe,NTU,T5);

end

state{3}.T = T3;
state{3}.s = janaf('s','N2',T3)*mN2_a + janaf('s','O2',T3)*mO2_a - s0a;

% States 4 and 5
state{5}.T = T5;

[mCO2_g,mH2O_g,mO2_g,mN2_g] = GasMassFraction(lambda,x,y,z);

h0g = mCO2_g*janaf('h','CO2',273.15) + mH2O_g*janaf('h','H2O',273.15)...
    + mO2_g*janaf('h','O2',273.15) + mN2_g*janaf('h','N2',273.15);

state{4}.h = mCO2_g*janaf('h','CO2',T4) + mH2O_g*janaf('h','H2O',T4)...
    + mO2_g*janaf('h','O2',T4) + mN2_g*janaf('h','N2',T4) - h0g;

state{5}.h = mCO2_g*janaf('h','CO2',T5) + mH2O_g*janaf('h','H2O',T5)...
    + mO2_g*janaf('h','O2',T5) + mN2_g*janaf('h','N2',T5) - h0g;


s0g = mCO2_g*janaf('s','CO2',273.15) + mH2O_g*janaf('s','H2O',273.15)...
    + mO2_g*janaf('s','O2',273.15) + mN2_g*janaf('s','N2',273.15);

state{4}.s = mCO2_g*janaf('s','CO2',T4) + mH2O_g*janaf('s','H2O',T4)...
    + mO2_g*janaf('s','O2',T4) + mN2_g*janaf('s','N2',T4) - s0g;

state{5}.s = mCO2_g*janaf('s','CO2',T5) + mH2O_g*janaf('s','H2O',T5)...
    + mO2_g*janaf('s','O2',T5) + mN2_g*janaf('s','N2',T5) - s0g;

%%%% State 6 %%%%

T = linspace(state{2}.T,state{3}.T,150);
cpa = mean(mO2_a*janaf('c','O2',T) + mN2_a*janaf('c','N2',T))*10^3;
T6_0 = T5 - 100;

% first iteration
T = linspace(T5,T6_0,150);
cpg = mean(mCO2_g*janaf('c','CO2',T) + mH2O_g*janaf('c','H2O',T)...
    + mO2_g*janaf('c','O2',T) + mN2_g*janaf('c','N2',T))*10^3;

T6 = T5 - ((ma*cpa)/(mg*cpg))*(state{3}.T-state{2}.T);

T6_prev = T6_0;

while abs(T6 - T6_prev) > 0.01
   
    T6_prev = T6;
    
    T = linspace(T5,T6,150);
    cpg = mean(mCO2_g*janaf('c','CO2',T) + mH2O_g*janaf('c','H2O',T)...
        + mO2_g*janaf('c','O2',T) + mN2_g*janaf('c','N2',T))*10^3;

    T6 = T5 - ((ma*cpa)/(mg*cpg))*(state{3}.T-state{2}.T);
end

state{6}.T = T6;

state{6}.h = mCO2_g*janaf('h','CO2',T6) + mH2O_g*janaf('h','H2O',T6)...
    + mO2_g*janaf('h','O2',T6) + mN2_g*janaf('h','N2',T6) - h0g;

state{6}.s = mCO2_g*janaf('s','CO2',T6) + mH2O_g*janaf('s','H2O',T6)...
    + mO2_g*janaf('s','O2',T6) + mN2_g*janaf('s','N2',T6) - s0g;


%%%% Exergy calculation %%%%
[state] = GasTurbineExergy(state);

%%%%%%%%%%%%%%%%%%
%%%% Analysis %%%%
%%%%%%%%%%%%%%%%%%

%%%% Energy analysis %%%%
[Energy_losses,labels_Energy,etaMec,etaCyclen,etaToten] = energyAnalysis6(state,ma,mc,mg,Pe,x,y,z);

%%%% Exergy analysis %%%%
[Exergy_losses,labels_Ex,etaMec,etaRotex,etaCyclex,etaCombex,etaTotex] = exergyAnalysis6(state,ma,mc,mg,x,y,z,Pe);

end