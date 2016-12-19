function [Exergy_losses,labels_Ex,etaMec,etaRotex,etaCyclex,etaCombex,etaTotex] = exergyAnalysis4(state,ma,mc,mg,x,y,z,Pe,fuel)
%Function calculating the exergy losses and exergy efficiencies of a gas
% turbine cycle without recuperator.
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
%   - Exergy_losses : a matrix with the exergy losses of the cycle in [MW]
%   - labels_Ex : a matrix with the labels associate at each loss of
%                     Exergy_losses
%   - etaMec : the mechanical efficiency of the cycle
%   - etaRotex : the efficiency of the turbo-machinery group
%   - etaCyclex : the cycle's energy efficiency
%   - etaCombex : the combustion efficiency
%   - etaTotex : the total energy efficiency of the cycle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% exergy of the fuel %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ec = fuel_exergy(x,y,z,fuel); % [kJ/kg]

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% work and power %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
WmT = state{3}.h - state{4}.h; % turbine work [kJ/kg]
WmC = state{2}.h - state{1}.h; % compressor work [kJ/kg]

PmT = mg*WmT; % turbine power [kW]
PmC = ma*WmC; % compressor power [kW]

Pm = PmT-PmC; % motor power [kW]
Pprim = mc*ec; % primary power [kW]

%%%%%%%%%%%%%%%%%%%%%%%
%%%% exergy losses %%%%
%%%%%%%%%%%%%%%%%%%%%%%

Pfmec = Pm-Pe; % mechanical losses power [kW]
irrComb = ec*mc + state{2}.e*ma - state{3}.e*mg; % irreversibilities at combustion in [kW]
irrTC = mg*(state{3}.e - state{4}.e) - ma*(state{2}.e - state{1}.e) - Pm; % irreversibilities at the turbo-machinery group [kW]
Pech = mg*state{4}.e; % Exhaust losses in [kW]

%%%%%%%%%%%%%%%%%%%%%%
%%%% efficiencies %%%%
%%%%%%%%%%%%%%%%%%%%%%

etaMec = Pe/Pm; % mechanical efficiency
etaRotex = Pm/(mg*(state{3}.e - state{4}.e) - ma*(state{2}.e - state{1}.e)); % efficiency of the turbo-machinery grou)
etaCyclex = Pm/(mg*state{3}.e - ma*state{2}.e); % efficiency of the cycle
etaCombex = (mg*state{3}.e - ma*state{2}.e)/(mc*ec); % efficiency of the combustion
etaTotex = Pe/Pprim; % total efficiency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% pie chart (in [MW]) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lab_EP = vertcat({'Effective power'},strcat({num2str(0.1*round(10*Pe/10^3))},{' '},{'MW'}));
lab_Mec = vertcat({'Mechanical losses'},strcat({num2str(0.1*round(10*Pfmec/10^3))},{' '},{'MW'}));
lab_Comb = vertcat({'Irreversibilities at combustion'},strcat({num2str(0.1*round(10*irrComb/10^3))},{' '},{'MW'}));
lab_TC = vertcat({'Turbo-machinery losses'},strcat({num2str(0.1*round(10*irrTC/10^3))},{' '},{'MW'}));
lab_Ex = vertcat({'Exhaust losses'},strcat({num2str(0.1*round(10*Pech/10^3))},{' '},{'MW'}));
labels_Ex = {lab_Mec,lab_Comb,lab_TC,lab_Ex,lab_EP};
Exergy_losses = [Pfmec/10^3 irrComb/10^3 irrTC/10^3 Pech/10^3 Pe/10^3];


end