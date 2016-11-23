function [Energy_losses,labels_Energy,etaMec,etaCyclen,etaToten] = energyAnalysis4(state,ma,mc,mg,Pe,x,y,z)
% INPUTS :
%   - state : cell array containing all the state characteristics with T in
%             [K], p in [bar], h in [kJ/kg], s in [kJ/(Kg*K)] and e in [kJ/kg]
%   - ma : air mass flow rate in [kg/s]
%   - mc : fuel mass flow rate in [kg/s]
%   - mg : gas mass flow rate in [kg/s]
%   - Pe : effective power in [kW]

%%%% Fuel characteristics %%%%
Mfuel = x*12+y+z*16; % molar mass
DrH = 282.4*((2*z+y*x-y)/(2+y)) + ((2*x-2*z+y)/(2+y))*((1-y/4)*393.4 + (y/4)*802.4);
LHV = (DrH/Mfuel)*10^3; %[kJ/kg] --> lower heating value

%%%% Energy analysis %%%%

WmT = state{3}.h - state{4}.h; % turbine work [kJ/kg]
WmC = state{2}.h - state{1}.h; % compressor work [kJ/kg]

PmT = mg*WmT; % turbine power [kW]
PmC = ma*WmC; % compressor power [kW]

Pm = PmT-PmC; % motor power [kW]
Pprim = mc*LHV; % primary power [kW]

Qcomb = (mc/ma)*LHV; % [kJ/kg_air]

Wm = Pm/ma; % motor work [kJ/kg_air]

Pfmec = Pm-Pe; % mechanical losses power [kW]
Pech = Pprim - Pe - Pfmec; % exhaust power [kW]

% efficiencies

etaMec = Pe/Pm;
etaCyclen = Wm/Qcomb;
etaToten = Pe/(mc*LHV);

% pie chart (in [MW])
Energy_losses = [Pe/10^3 Pech/10^3 Pfmec/10^3];
labels_Energy = {'Effective power','Exhaust losses','Mechanical losses'};

end