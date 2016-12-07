function [Exergy_losses,labels_Ex,etaMec,etaRotex,etaCyclex,etaCombex,etaTotex] = exergyAnalysis4(state,ma,mc,mg,x,y,z,Pe,fuel)


%%%% exergy of the fuel %%%%
ec = fuel_exergy(x,y,z,fuel); % [kJ/kg]

%%%% work and power %%%%
WmT = state{3}.h - state{4}.h; % turbine work [kJ/kg]
WmC = state{2}.h - state{1}.h; % compressor work [kJ/kg]

PmT = mg*WmT; % turbine power [kW]
PmC = ma*WmC; % compressor power [kW]

Pm = PmT-PmC; % motor power [kW]
Pprim = mc*ec; % primary power [kW]

%%%% exergy losses %%%%

Pfmec = Pm-Pe; % mechanical losses power [kW]
irrComb = ec*mc + state{2}.e*ma - state{3}.e*mg;
irrTC = mg*(state{3}.e - state{4}.e) - ma*(state{2}.e - state{1}.e) - Pm;
Pech = mg*state{4}.e;

%%%% efficiencies %%%%

etaMec = Pe/Pm;
etaRotex = Pm/(mg*(state{3}.e - state{4}.e) - ma*(state{2}.e - state{1}.e));
etaCyclex = Pm/(mg*state{3}.e - ma*state{2}.e);
etaCombex = (mg*state{3}.e - ma*state{2}.e)/(mc*ec);
etaTotex = Pe/Pprim;

%%%% pie chart (in [MW]) %%%%
lab_EP = vertcat({'Effective power'},strcat({num2str(0.1*round(10*Pe/10^3))},{' '},{'MW'}));
lab_Mec = vertcat({'Mechanical losses'},strcat({num2str(0.1*round(10*Pfmec/10^3))},{' '},{'MW'}));
lab_Comb = vertcat({'Irreversibilities at combustion'},strcat({num2str(0.1*round(10*irrComb/10^3))},{' '},{'MW'}));
lab_TC = vertcat({'Turbo-machinery losses'},strcat({num2str(0.1*round(10*irrTC/10^3))},{' '},{'MW'}));
lab_Ex = vertcat({'Exhaust losses'},strcat({num2str(0.1*round(10*Pech/10^3))},{' '},{'MW'}));
labels_Ex = {lab_Mec,lab_Comb,lab_TC,lab_Ex,lab_EP};
Exergy_losses = [Pfmec/10^3 irrComb/10^3 irrTC/10^3 Pech/10^3 Pe/10^3];


end