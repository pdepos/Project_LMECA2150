function [Energy_losses,labels_Energy,etaMec,etaCyclen,etaToten] = energyAnalysis6(state,ma,mc,mg,Pe,x,y,z,fuel)
% Function calculating the energy losses and energy efficiencies of a gas
% turbine cycle with a recuperator.
%
% INPUTS :
%   - state : cell array containing all the state characteristics with T in
%             [K], p in [bar], h in [kJ/kg], s in [kJ/(Kg*K)] and e in [kJ/kg]
%   - ma : air mass flow rate in [kg/s]
%   - mc : fuel mass flow rate in [kg/s]
%   - mg : gas mass flow rate in [kg/s]
%   - Pe : effective power in [kW]
%   - x : the number of carbon of the fuel molecule (CxHyOz)
%   - y : the number of hydrogen of the fuel molecule (CxHyOz)
%   - z : the number of oxygen of the fuel molecule (CxHyOz)
%   - fuel : the fuel formula in string form. It MUST be 'CH4' or 'C12H23'
%
% OUTPUTS :
%   - Energy_losses : a matrix with the energy losses of the cycle in [MW]
%   - labels_Energy : a matrix with the labels associate at each loss of
%                     Energy_losses
%   - etaMec : the mechanical efficiency of the cycle
%   - etaCyclen : the cycle's energy efficiency
%   - etaToten : the total energy efficiency of the cycle


%%%% Fuel characteristics %%%%

if strcmp(fuel,'CH4')
    LHV = 50.1*10^3; %[kJ/kg]
elseif strcmp(fuel,'C12H23')
    LHV = 41.76*10^3; %[kJ/kg]
end

%%%% Energy analysis %%%%
WmT = state{4}.h - state{5}.h; % turbine work [kJ/kg]
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
lab_EP = vertcat({'Effective power'},strcat({num2str(0.1*round(10*Pe/10^3))},{' '},{'MW'}));
lab_ExL = vertcat({'Exhaust losses'},strcat({num2str(0.1*round(10*Pech/10^3))},{' '},{'MW'}));
lab_Mec = vertcat({'Mechanical losses'},strcat({num2str(0.1*round(10*Pfmec/10^3))},{' '},{'MW'}));
labels_Energy = {lab_EP,lab_ExL,lab_Mec};



end