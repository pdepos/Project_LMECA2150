function [stateV,stateTG,Energy_lossesTGV,labels_EnergyTGV,Exergy_lossesTGV,labels_ExTGV,massFlow,energyEff,exergyEff] = ...
    mainCombined3(Pe_gas,pHP,pMP,pLP,Twater,Dta,Dtp_LP,Dtp_MP,Dtp_HP,x,y,z,T1g,r,kcc,T3g,fuel)

%%%%%%%%%%%%%%
%%%% DATA %%%%
%%%%%%%%%%%%%%

%%%% isentropic efficiencies %%%%

etaSiP1 = 0.9; 
etaSiP2 = 0.9;
etaSiP3 = 0.9;
etaSiHP = 0.88;
etaSiMP = 0.88;
etaSiLP = 0.88;

%%%% polytropic efficiencies (Gas turbine) %%%%
etaPiT = 0.9;
etaPiC = 0.9;

%%%% other parameters %%%%
Tpinch_water = 15;
x20 = 0.9;

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
stateV = State_creation_Combine(20);

%%%% State 1 (satured liquid) %%%%
stateV{1}.T = Twater + Tpinch_water;
stateV{1}.p = XSteam('psat_T',stateV{1}.T);
stateV{1}.x = 0;
stateV{1}.h = XSteam('hL_T',stateV{1}.T);
stateV{1}.s = XSteam('sL_T',stateV{1}.T);

%%%% State 14 (overheated steam) %%%%
stateV{14}.p = pHP;
stateV{14}.T = stateTG{4}.T - Dta;
stateV{14}.h = XSteam('h_pT',stateV{14}.p,stateV{14}.T);
stateV{14}.s = XSteam('s_pT',stateV{14}.p,stateV{14}.T);

%%%% State 11 (satured steam) %%%%
stateV{11}.p = stateV{14}.p;
stateV{11}.T = XSteam('Tsat_p',stateV{11}.p);
stateV{11}.h = XSteam('hV_T',stateV{11}.T);
stateV{11}.s = XSteam('sV_T',stateV{11}.T);
stateV{11}.x = 1;

%%%% State 10 (satured liquid) %%%%
stateV{10}.T = stateV{11}.T;
stateV{10}.p = stateV{11}.p;
stateV{10}.x = 0;
stateV{10}.h = XSteam('hL_T',stateV{10}.T);
stateV{10}.s = XSteam('sL_T',stateV{10}.T);

%%%% isobaric valves %%%%
stateV{17}.p = pMP;
stateV{16}.p = stateV{17}.p;
stateV{15}.p = stateV{16}.p;
stateV{13}.p = stateV{16}.p;

stateV{19}.p = pLP;
stateV{18}.p = stateV{19}.p;
stateV{6}.p = stateV{19}.p;

%%%% isobaric exchanger %%%%
stateV{5}.p = stateV{6}.p;
stateV{4}.p = stateV{5}.p;
stateV{3}.p = stateV{4}.p;
stateV{2}.p = stateV{3}.p;

%%%% State 15 (overheated steam) %%%%
s15s = stateV{14}.s;
h15s = XSteam('h_ps',stateV{15}.p,s15s);

stateV{15}.h = stateV{14}.h - etaSiHP*(stateV{14}.h - h15s);
stateV{15}.T = XSteam('T_ph',stateV{15}.p,stateV{15}.h);
stateV{15}.s = XSteam('s_pT',stateV{15}.p,stateV{15}.T);

%%%% State 12 (satured steam) %%%%
stateV{12}.p = stateV{13}.p;
stateV{12}.T = XSteam('Tsat_p',stateV{12}.p);
stateV{12}.h = XSteam('hV_T',stateV{12}.T);
stateV{12}.s = XSteam('sV_T',stateV{12}.T);
stateV{12}.x = 1;

%%%% State 8 (satured liquid) %%%%
stateV{8}.p = stateV{12}.p;
stateV{8}.T = XSteam('Tsat_p',stateV{8}.p);
stateV{8}.h = XSteam('hL_T',stateV{8}.T);
stateV{8}.s = XSteam('sL_T',stateV{8}.T);
stateV{8}.x = 0;

%%%% State 9 (subcooled liquid) %%%%
stateV{9}.p = stateV{10}.p;
s9s = stateV{8}.s;
h9s = XSteam('h_ps',stateV{9}.p,s9s);

stateV{9}.h = stateV{8}.h + etaSiP3*(h9s - stateV{8}.h);
stateV{9}.T = XSteam('T_ph',stateV{9}.p,stateV{9}.h);
stateV{9}.s = XSteam('s_pT',stateV{9}.p,stateV{9}.T);

%%%% State 4 (satured steam) %%%%
stateV{4}.T = XSteam('Tsat_p',stateV{4}.p);
stateV{4}.h = XSteam('hV_T',stateV{4}.T);
stateV{4}.s = XSteam('sV_T',stateV{4}.T);
stateV{4}.x = 1;

%%%% State 3 (satured liquid) %%%%
stateV{3}.T = XSteam('Tsat_p',stateV{3}.p);
stateV{3}.h = XSteam('hL_T',stateV{3}.T);
stateV{3}.s = XSteam('sL_T',stateV{3}.T);
stateV{3}.x = 0;

%%%% State 2 (subcooled liquid) %%%%
s2s = stateV{1}.s;
h2s = XSteam('h_ps',stateV{2}.p,s2s);

stateV{2}.h = stateV{1}.h + etaSiP1*(h2s - stateV{1}.h);
stateV{2}.T = XSteam('T_ph',stateV{2}.p,stateV{2}.h);
stateV{2}.s = XSteam('s_pT',stateV{2}.p,stateV{2}.T);

%%%% State 7 (subcooled liquid) %%%%
stateV{7}.p = stateV{8}.p;
s7s = stateV{3}.s;
h7s = XSteam('h_ps',stateV{7}.p,s7s);

stateV{7}.h = stateV{3}.h + etaSiP2*(h7s - stateV{3}.h);
stateV{7}.T = XSteam('T_ph',stateV{7}.p,stateV{7}.h);
stateV{7}.s = XSteam('s_pT',stateV{7}.p,stateV{7}.T);

%%%% State 17 (overheated steam) %%%%
stateV{17}.T = stateTG{4}.T - Dta;
stateV{17}.h = XSteam('h_pT',stateV{17}.p,stateV{17}.T);
stateV{17}.s = XSteam('s_pT',stateV{17}.p,stateV{17}.T);

%%%% State 18 (overheated steam) %%%%
s18s = stateV{17}.s;
h18s = XSteam('h_ps',stateV{18}.p,s18s);

stateV{18}.h = stateV{17}.h - etaSiMP*(stateV{17}.h - h18s);
stateV{18}.T = XSteam('T_ph',stateV{18}.p,stateV{18}.h);
stateV{18}.s = XSteam('s_pT',stateV{18}.p,stateV{18}.T);

%%%% State 13 (overheated steam) %%%%
stateV{13}.T = stateV{10}.T;
stateV{13}.h = XSteam('h_pT',stateV{13}.p,stateV{13}.T);
stateV{13}.s = XSteam('s_pT',stateV{13}.p,stateV{13}.T);

%%%% State 5 (overheated steam) %%%%
stateV{5}.T = stateV{8}.T;
stateV{5}.h = XSteam('h_pT',stateV{5}.p,stateV{5}.T);
stateV{5}.s = XSteam('s_pT',stateV{5}.p,stateV{5}.T);


%%%% State 6 (overheated steam) %%%%
stateV{6}.T = stateV{10}.T;
stateV{6}.h = XSteam('h_pT',stateV{6}.p,stateV{6}.T);
stateV{6}.s = XSteam('s_pT',stateV{6}.p,stateV{6}.T);

%%%% Pinch points %%%%

TpHP = stateV{11}.T + Dtp_HP;
TpMP = stateV{5}.T + Dtp_MP;
TpLP = stateV{4}.T + Dtp_LP;

[mCO2_g,mH2O_g,mO2_g,mN2_g,Mg] = GasMassFraction(lambda,x,y,z);

h0g = mCO2_g*janaf('h','CO2',273.15) + mH2O_g*janaf('h','H2O',273.15)...
    + mO2_g*janaf('h','O2',273.15) + mN2_g*janaf('h','N2',273.15);

hpLP = mCO2_g*janaf('h','CO2',TpLP+273.15) + mH2O_g*janaf('h','H2O',TpLP+273.15)...
    + mO2_g*janaf('h','O2',TpLP+273.15) + mN2_g*janaf('h','N2',TpLP+273.15) - h0g;
hpMP = mCO2_g*janaf('h','CO2',TpMP+273.15) + mH2O_g*janaf('h','H2O',TpMP+273.15)...
    + mO2_g*janaf('h','O2',TpMP+273.15) + mN2_g*janaf('h','N2',TpMP+273.15) - h0g;
hpHP = mCO2_g*janaf('h','CO2',TpHP+273.15) + mH2O_g*janaf('h','H2O',TpHP+273.15)...
    + mO2_g*janaf('h','O2',TpHP+273.15) + mN2_g*janaf('h','N2',TpHP+273.15) - h0g;

%%%% Mass flows and mixing valves %%%%

% the unknown vector is x = [mvHP mvMP mvLP h16]

mf = @(x)([mg*(stateTG{4}.h-hpHP)-x(1)*(stateV{14}.h-stateV{10}.h)-(x(1)+x(2))*(stateV{17}.h-x(4));...
    mg*(hpHP-hpMP)-x(3)*(stateV{6}.h-stateV{5}.h) - x(2)*(stateV{13}.h-stateV{8}.h) - x(1)*(stateV{10}.h-stateV{9}.h);...
    mg*(hpMP-hpLP)-(x(1)+x(2))*(stateV{8}.h-stateV{7}.h) - x(3)*(stateV{5}.h-stateV{3}.h);...
    (x(1)+x(2))*x(4) - x(1)*stateV{15}.h - x(2)*stateV{13}.h]);

opts = optimoptions(@fsolve,'Display','none');
Soluce = fsolve(mf,[80 11 10 400],opts);

mvHP = Soluce(1);
mvMP = Soluce(2);
mvLP = Soluce(3);
h16 = Soluce(4);

%%%% State 16 (overheated steam) %%%%
stateV{16}.h = h16;
stateV{16}.T = XSteam('T_ph',stateV{16}.p,stateV{16}.h);
stateV{16}.s = XSteam('s_pT',stateV{16}.p,stateV{16}.T);

%%%% State 19 (overheated steam) %%%%
stateV{19}.p = pLP;
stateV{19}.h = ((mvHP+mvMP)*stateV{18}.h + mvLP*stateV{6}.h)/(mvHP+mvMP+mvLP);
stateV{19}.T = XSteam('T_ph',stateV{19}.p,stateV{19}.h);
stateV{19}.s = XSteam('s_ph',stateV{19}.p,stateV{19}.h);

%%%% State 20 (satured steam) %%%%
stateV{20}.T = Twater + Tpinch_water;
stateV{20}.p = XSteam('Tsat_p',stateV{20}.T);

s20s = stateV{19}.s;
x20s = (s20s - XSteam('sL_T',stateV{20}.T))/(XSteam('sV_T',stateV{20}.T) - XSteam('sL_T',stateV{20}.T));
h20s = x20s*XSteam('hV_T',stateV{20}.T) + (1-x20s)*XSteam('hL_T',stateV{20}.T);

stateV{20}.h = stateV{19}.h - etaSiLP*(stateV{19}.h - h20s);
stateV{20}.x = (stateV{20}.h - XSteam('hL_T',stateV{20}.T))/(XSteam('hV_T',stateV{20}.T) - XSteam('hL_T',stateV{20}.T));
stateV{20}.s = stateV{20}.x*XSteam('sV_T',stateV{20}.T) + (1-stateV{20}.x)*XSteam('sL_T',stateV{20}.T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Exhaust point (Gas cycle) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h5g = hpLP - ((mvLP+mvMP+mvHP)*(stateV{3}.h-stateV{2}.h))/mg;
stateTG{5}.h = h5g;

T5g = Temperature5g(lambda,h5g,x,y,z);
stateTG{5}.T = T5g;

s0g = mCO2_g*janaf('s','CO2',273.15) + mH2O_g*janaf('s','H2O',273.15)...
    + mO2_g*janaf('s','O2',273.15) + mN2_g*janaf('s','N2',273.15);

stateTG{5}.s = mCO2_g*janaf('s','CO2',T5g+273.15) + mH2O_g*janaf('s','H2O',T5g+273.15)...
    + mO2_g*janaf('s','O2',T5g+273.15) + mN2_g*janaf('s','N2',T5g+273.15) - s0g;

stateTG(5) = GasTurbineExergy(stateTG(5));
stateTG{5}.p = stateTG{4}.p;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Exergy calculation (steam cycle) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stateV] = VapExergy(stateV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Energy efficiencies %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Steam cycle %%%%

[Pcond,PfmecV,etaMecV,etaCyclenV,etaTotenV] = Vap_energyAnalysis3(stateV,mvHP,mvMP,mvLP);

% motor power
PmT = mvHP*(stateV{14}.h-stateV{15}.h) + (mvHP+mvMP)*(stateV{17}.h-stateV{18}.h)...
    + (mvHP+mvMP+mvLP)*(stateV{19}.h-stateV{20}.h); % [kW]

%%%% Electrical power %%%%
PeV = PmT - PfmecV; %[kW]

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

%%%% Steam Cycle %%%%

% Condenser losses
Pcondex = (mvHP+mvMP+mvLP)*(stateV{20}.e - stateV{1}.e); %[kW]

% Turbines losses
irrT = mvHP*(stateV{14}.e-stateV{15}.e) + (mvHP+mvMP)*(stateV{17}.e - stateV{18}.e)...
    + (mvHP+mvMP+mvLP)*(stateV{19}.e-stateV{20}.e) - PmT; % [kW]

%%%% Heat Exchanger %%%%
irrHE = abs(mg*(stateTG{5}.e - stateTG{4}.e) + mvHP*(stateV{14}.e-stateV{9}.e)...
    + mvMP*(stateV{13}.e-stateV{8}.e) + mvLP*(stateV{6}.e-stateV{3}.e)...
    + (mvHP+mvMP+mvLP)*(stateV{3}.e-stateV{2}.e)...
    + (mvHP+mvMP)*(stateV{17}.e-stateV{16}.e+stateV{8}.e-stateV{7}.e)); %[kW]

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
etaCyclexV = PmT/( mvHP*(stateV{14}.e-stateV{9}.e)+ mvMP*(stateV{13}.e-stateV{8}.e)...
    + mvLP*(stateV{6}.e-stateV{3}.e) + (mvHP+mvMP+mvLP)*(stateV{3}.e-stateV{2}.e)...
    + (mvHP+mvMP)*(stateV{17}.e-stateV{16}.e+stateV{8}.e-stateV{7}.e));
etaTransex = ( mvHP*(stateV{14}.e-stateV{9}.e) + mvMP*(stateV{13}.e-stateV{8}.e)...
    + mvLP*(stateV{6}.e-stateV{3}.e) + (mvHP+mvMP+mvLP)*(stateV{3}.e-stateV{2}.e)...
    + (mvHP+mvMP)*(stateV{17}.e-stateV{16}.e+stateV{8}.e-stateV{7}.e))/(mg*(stateTG{4}.e - stateTG{5}.e));
etaTotexV = PeV/(mg*(stateTG{4}.e - stateTG{5}.e));

%%%% Combined exergy efficiency %%%%
ec = fuel_exergy(x,y,z,fuel);
etaTotexTGV = (PeTGV*10^3)/(ec*mc);

%%%%%%%%%%%%%%%%%
%%%% Results %%%%
%%%%%%%%%%%%%%%%%

massFlow = [ma,mc,mg,mvHP,mvMP,mvLP];
energyEff = [etaMecV,etaCyclenV,etaTotenV,etaTotenTGV,etaMecTG,etaCyclenTG,etaTotenTG];
exergyEff = [etaMecTG,etaRotexTG,etaCyclexTG,etaCombexTG,etaTotexTG,...
    etaRotexV,etaCyclexV,etaTransex,etaTotexV,etaTotexTGV];
end