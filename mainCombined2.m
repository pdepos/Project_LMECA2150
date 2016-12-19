function [stateV,stateTG,Energy_lossesTGV,labels_EnergyTGV,Exergy_lossesTGV,labels_ExTGV,massFlow,energyEff,exergyEff] = mainCombined2(Pe_gas,pHP,pLP,Twater,Dta,Dtp_LP,Dtp_HP,x,y,z,T1g,r,kcc,T3g,fuel)



%%%%%%%%%%%%%%
%%%% DATA %%%%
%%%%%%%%%%%%%%

%%%% isentropic efficiencies %%%%

etaSiP1 = 0.9; 
etaSiP2 = 0.9;
etaSiHP = 0.9;
etaSiLP = 0.9;

%%%% polytropic efficiencies (Gas turbine) %%%%
etaPiT = 0.9;
etaPiC = 0.9;

%%%% other parameters %%%%
Tpinch_water = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gas turbine calculation %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stateG,Energy_lossesTG,labels_EnergyTG,etaMecTG,etaCyclenTG,etaTotenTG,Exergy_lossesTG,labels_ExTG,etaRotexTG,etaCyclexTG,etaCombexTG,etaTotexTG,ma,mc,mg,lambda] = mainTurbineGaz(T1g,r,etaPiC,kcc,T3g,etaPiT,Pe_gas,x,y,z,fuel);

% adding the fift state of the gas (exhaust)
stateTG = State_creation_Turbine(5);
for i = 1:4
   stateG{i}.T = stateG{i}.T - 273.15;
   stateTG{i} = stateG{i}; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% States calculation (vapor cycle) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% State creation %%%%
stateV = State_creation_Combine(10);

%%%% State 1 (satured liquid) %%%%
stateV{1}.T = Twater + Tpinch_water;
stateV{1}.p = XSteam('psat_T',stateV{1}.T);
stateV{1}.x = 0;
stateV{1}.h = XSteam('hL_T',stateV{1}.T);
stateV{1}.s = XSteam('sL_T',stateV{1}.T);

%%%% State 8 (overheating vapor) %%%%
stateV{8}.T = stateTG{4}.T - Dta;
stateV{8}.p = pHP;
stateV{8}.h = XSteam('h_pT',stateV{8}.p,stateV{8}.T);
stateV{8}.s = XSteam('s_pT',stateV{8}.p,stateV{8}.T);

%%%% State 7 (satured vapor) %%%%
stateV{7}.p = stateV{8}.p;
stateV{7}.T = XSteam('Tsat_p',stateV{7}.p);
stateV{7}.h = XSteam('hV_T',stateV{7}.T);
stateV{7}.s = XSteam('sV_T',stateV{7}.T);
stateV{7}.x = 1;

%%%% State 6 (satured liquid) %%%%
stateV{6}.p = stateV{8}.p;
stateV{6}.T = stateV{7}.T;
stateV{6}.x = 0;
stateV{6}.h = XSteam('hL_T',stateV{6}.T);
stateV{6}.s = XSteam('sL_T',stateV{6}.T);

%%%% State 9 (overheating vapor) %%%%
stateV{9}.p = pLP;

%%%% State 2 (subcooled liquid) %%%%
stateV{2}.p = stateV{9}.p;

s2s = stateV{1}.s;
h2s = XSteam('h_ps',stateV{2}.p,s2s);

stateV{2}.h = stateV{1}.h + etaSiP1*(h2s - stateV{1}.h);
stateV{2}.T = XSteam('T_ph',stateV{2}.p,stateV{2}.h);
stateV{2}.s = XSteam('s_pT',stateV{2}.p,stateV{2}.T);

%%%% State 3 (satured liquid) %%%%
stateV{3}.p = stateV{9}.p;
stateV{3}.T = XSteam('Tsat_p',stateV{3}.p);
stateV{3}.x = 0;
stateV{3}.h = XSteam('hL_T',stateV{3}.T);
stateV{3}.s = XSteam('sL_T',stateV{3}.T);

%%%% State 4 (satured vapor) %%%%
stateV{4}.p = stateV{9}.p;
stateV{4}.T = XSteam('Tsat_p',stateV{4}.p);
stateV{4}.x = 1;
stateV{4}.h = XSteam('hV_T',stateV{4}.T);
stateV{4}.s = XSteam('sV_T',stateV{4}.T);

%%%% State 5 (subcooled liquid) %%%%
stateV{5}.p = stateV{6}.p;

s5s = stateV{3}.s;
h5s = XSteam('h_ps',stateV{5}.p,s5s);

stateV{5}.h = stateV{3}.h + etaSiP2*(h5s - stateV{3}.h);
stateV{5}.T = XSteam('T_ph',stateV{5}.p,stateV{5}.h);
stateV{5}.s = XSteam('s_pT',stateV{5}.p,stateV{5}.T);

%%%% State 9 (overheating vapor) %%%%
s9s = stateV{8}.s;
h9s = XSteam('h_ps',stateV{9}.p,s9s);

stateV{9}.h = stateV{8}.h - etaSiHP*(stateV{8}.h - h9s);
stateV{9}.T = XSteam('T_ph',stateV{9}.p,stateV{9}.h);
stateV{9}.s = XSteam('s_pT',stateV{9}.p,stateV{9}.T);

%%%% State 10 (satured vapor) %%%%
stateV{10}.T = Twater + Tpinch_water;
stateV{10}.p = XSteam('psat_T',stateV{10}.T);

s10s = stateV{9}.s;
x10s = (s10s - XSteam('sL_T',stateV{10}.T))/(XSteam('sV_T',stateV{10}.T) - XSteam('sL_T',stateV{10}.T));
h10s = x10s*XSteam('hV_T',stateV{10}.T) + (1-x10s)*XSteam('hL_T',stateV{10}.T);

stateV{10}.h = stateV{9}.h - etaSiLP*(stateV{9}.h - h10s);
stateV{10}.x = (stateV{10}.h - XSteam('hL_T',stateV{10}.T))/(XSteam('hV_T',stateV{10}.T) - XSteam('hL_T',stateV{10}.T));
stateV{10}.s = stateV{10}.x*XSteam('sV_T',stateV{10}.T) + (1-stateV{10}.x)*XSteam('sL_T',stateV{10}.T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Mass flow rates calculation (vapor cycle) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Pinch points %%%%

% LP pinch point
TpLP = stateV{3}.T + Dtp_LP;

[mCO2_g,mH2O_g,mO2_g,mN2_g,Mg] = GasMassFraction(lambda,x,y,z);

h0g = mCO2_g*janaf('h','CO2',273.15) + mH2O_g*janaf('h','H2O',273.15)...
    + mO2_g*janaf('h','O2',273.15) + mN2_g*janaf('h','N2',273.15);

hpLP = mCO2_g*janaf('h','CO2',TpLP+273.15) + mH2O_g*janaf('h','H2O',TpLP+273.15)...
    + mO2_g*janaf('h','O2',TpLP+273.15) + mN2_g*janaf('h','N2',TpLP+273.15) - h0g;

% HP pinch point
TpHP = stateV{6}.T + Dtp_HP;

hpHP = mCO2_g*janaf('h','CO2',TpHP+273.15) + mH2O_g*janaf('h','H2O',TpHP+273.15)...
    + mO2_g*janaf('h','O2',TpHP+273.15) + mN2_g*janaf('h','N2',TpHP+273.15) - h0g;

%%%% mass flow rates %%%%

mvHP = (mg*(stateTG{4}.h-hpHP))/(stateV{8}.h-stateV{6}.h);
mvLP = (mg*(hpHP-hpLP)-mvHP*(stateV{6}.h-stateV{5}.h))/(stateV{9}.h-stateV{3}.h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Exhaust point (Gas cycle) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h5g = hpLP - ((mvLP+mvHP)*(stateV{3}.h-stateV{2}.h))/mg;
stateTG{5}.h = h5g;
T5g = Temperature5g(lambda,h5g,x,y,z);
stateTG{5}.T = T5g;

s0g = mCO2_g*janaf('s','CO2',273.15) + mH2O_g*janaf('s','H2O',273.15)...
    + mO2_g*janaf('s','O2',273.15) + mN2_g*janaf('s','N2',273.15);

stateTG{5}.s = mCO2_g*janaf('s','CO2',T5g+273.15) + mH2O_g*janaf('s','H2O',T5g+273.15)...
    + mO2_g*janaf('s','O2',T5g+273.15) + mN2_g*janaf('s','N2',T5g+273.15) - s0g;

stateTG(5) = GasTurbineExergy(stateTG(5));
%stateTG{5}.e = (stateTG{5}.h-h0g)-(15+273.15)*(stateTG{5}.s-s0g);
stateTG{5}.p = stateTG{4}.p;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Exergy calculation (steam cycle) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stateV] = VapExergy(stateV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Energy efficiencies %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Gas Turbine %%%%
%[Energy_lossesTG,labels_EnergyTG,etaMecTG,etaCyclenTG,etaTotenTG] = energyAnalysis4(stateTG,ma,mc,mg,Pe_gas*10^3,x,y,z,fuel);
%%%% Steam cycle %%%%
[Pcond,PfmecV,etaMecV,etaCyclenV,etaTotenV] = Vap_energyAnalysis2(stateV,mvHP,mvLP); % losses given in [kW]

%%%% Electrical power %%%%
PeV = mvHP*(stateV{8}.h-stateV{9}.h) + (mvHP+mvLP)*(stateV{9}.h-stateV{10}.h) - PfmecV; %[kW]

%%%% Assembly %%%%
PfmecTGV = PfmecV/10^3 + Energy_lossesTG(3); % [MW]
PeTGV = Pe_gas + PeV/10^3; % [MW]

% losses at the stag
Pech = stateTG{5}.h*mg; %[kW]

%%%% Pie chart %%%%
Energy_lossesTGV = [PeTGV Pech/10^3 PfmecTGV Pcond/10^3];
lab_EP = vertcat({'Effective power'},strcat({num2str(0.1*round(10*PeTGV))},{' '},{'MW'}));
lab_ExL = vertcat({'Exhaust losses'},strcat({num2str(0.1*round(10*Pech/10^3))},{' '},{'MW'}));
lab_Mec = vertcat({'Mechanical losses'},strcat({num2str(0.1*round(10*PfmecTGV))},{' '},{'MW'}));
lab_cond = vertcat({'Losses at condenser'},strcat({num2str(0.1*round(10*Pcond/10^3))},{' '},{'MW'}));
labels_EnergyTGV = {lab_EP,lab_ExL,lab_Mec,lab_cond};

if strcmp(fuel,'CH4')
    LHV = 50.1*10^3; %[kJ/kg]
elseif strcmp(fuel,'C12H23')
    LHV = 41.76*10^3; %[kJ/kg]
end

etaTotenTGV = (PeTGV*10^3)/(LHV*mc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Exergy efficiencies %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PmT = mvHP*(stateV{8}.h-stateV{9}.h) + (mvHP+mvLP)*(stateV{9}.h-stateV{10}.h); % [kW]

%%%% Gas Turbine %%%%
%[Exergy_lossesTG,labels_ExTG,etaMecTG,etaRotexTG,etaCyclexTG,etaCombexTG,etaTotexTG] = exergyAnalysis4(stateTG,ma,mc,mg,x,y,z,Pe_gas*10^3,fuel);

%%%% Steam Cycle %%%%

% Condenser losses
Pcondex = (mvHP+mvLP)*(stateV{10}.e - stateV{1}.e); %[kW]

% Turbines losses
irrT = mvHP*(stateV{8}.e-stateV{9}.e) + (mvHP+mvLP)*(stateV{9}.e - stateV{10}.e) - PmT; % [kW]

%%%% Heat Exchanger %%%%
irrHE = abs(mg*(stateTG{5}.e - stateTG{4}.e) + mvHP*(stateV{8}.e-stateV{5}.e) +...
    mvLP*(stateV{9}.e-stateV{3}.e) + (mvHP+mvLP)*(stateV{3}.e-stateV{2}.e)); % [kW]

%%%% Stag losses %%%%
Pechex = mg*stateTG{5}.e; %[kW]

%%%% Combustion %%%%
irrComb = Exergy_lossesTG(2); %[MW]

%%%% Assembly %%%%
irrTC = irrT/10^3 + Exergy_lossesTG(3); %[MW]

%%%% pie chart (in [MW]) %%%%
lab_EP = vertcat({'Effective power'},strcat({num2str(0.1*round(10*PeTGV))},{' '},{'MW'}));
lab_Mec = vertcat({'Mechanical losses'},strcat({num2str(0.1*round(10*PfmecTGV))},{' '},{'MW'}));
lab_Comb = vertcat({'Irreversibilities'},{'at combustion'},strcat({num2str(0.1*round(10*irrComb))},{' '},{'MW'}));
lab_TC = vertcat({'Turbo-machinery losses'},strcat({num2str(0.1*round(10*irrTC))},{' '},{'MW'}));
lab_Ex = vertcat({'Exhaust losses'},strcat({num2str(0.1*round(10*Pechex/10^3))},{' '},{'MW'}));
lab_cond = vertcat({'Losses at condenser'},strcat({num2str(0.1*round(10*Pcondex/10^3))},{' '},{'MW'}));
lab_HE = vertcat({'Irreversibilities at'},{'the Heat Exchanger'},strcat({num2str(0.1*round(10*irrHE/10^3))},{' '},{'MW'}));
labels_ExTGV = {lab_Mec,lab_Comb,lab_TC,lab_Ex,lab_EP,lab_cond,lab_HE};
Exergy_lossesTGV = [PfmecTGV irrComb irrTC Pechex/10^3 PeTGV Pcondex/10^3 irrHE/10^3];

%%%% Exergy efficiencies (Steam cycle) %%%%

etaRotexV = PmT/(irrT + PmT);
etaCyclexV = PmT/(mvHP*(stateV{8}.e-stateV{5}.e) + mvLP*(stateV{9}.e-stateV{3}.e)...
    + (mvHP+mvLP)*(stateV{3}.e-stateV{2}.e));
etaTransex = (mvHP*(stateV{8}.e-stateV{5}.e) + mvLP*(stateV{9}.e-stateV{3}.e)...
    + (mvHP+mvLP)*(stateV{3}.e-stateV{2}.e))/(mg*(stateTG{4}.e - stateTG{5}.e));
etaTotexV = PeV/(mg*(stateTG{4}.e - stateTG{5}.e));

%%%% Combined exergy efficiency %%%%
ec = fuel_exergy(x,y,z,fuel);
etaTotexTGV = (PeTGV*10^3)/(ec*mc);

%%%%%%%%%%%%%%%%%
%%%% Results %%%%
%%%%%%%%%%%%%%%%%

massFlow = [ma,mc,mg,mvHP,mvLP];
energyEff = [etaMecV,etaCyclenV,etaTotenV,etaTotenTGV,etaMecTG,etaCyclenTG,etaTotenTG];
exergyEff = [etaMecTG,etaRotexTG,etaCyclexTG,etaCombexTG,etaTotexTG,...
    etaRotexV,etaCyclexV,etaTransex,etaTotexV,etaTotexTGV];
end