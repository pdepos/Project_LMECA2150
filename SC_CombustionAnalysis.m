function [ output_args ] = SC_CombustionAnalysis( x, y ,z, t_exh, lambda )
%UNTITLED Summary of this function goes here
%   Input Args:
%       - CxHyOz: Fuel type
%       - t_exh: Temperature of the exhaust gasses
%       - lambda: Excess Air Coefficient, cannot be smaller than 1

x = 1;
y = 4;
z = 0;
Ta = 15; %Ambiant Air Temperature
t_exh = 120;
lambda = 1.05;
LHVCH4 = 50150; %kj/kg_CH4

if lambda < 1
    lambda = 1;
end

MWO2 = 31.998; %kg/kmol
MWCO2 = 44.008;
MWH2O = 18.01494;
MWN2 = 28.014;

ma1 = ((MWO2+3.76*MWN2)*(x + (y - 2*z)/4))/(12*x + y + 16*z); % pouvoir comburivore [kg_air/kg_fuel]

Wfgasses      = (x*MWCO2 + 0.5*y*MWH2O + (lambda - 1)*MWO2 + (x + 0.25*(y-2*z))*3.76*MWN2); %kg
kmolesfgasses = (x + 0.5*y + (lambda - 1) + (x + 0.25*(y-2*z))*3.76); %kmoles
MWfgasses     = Wfgasses / kmolesfgasses; %kg_fg /kmol_fg
            
        
CO2_frac = x*MWCO2/Wfgasses; %kg_CO2 / kg_fg
H2O_frac = 0.5*y*MWH2O/Wfgasses; %kg_H2O / kg_fg
O2_frac  = (lambda - 1)*MWO2/Wfgasses;
N2_frac  = (x + 0.25*(y-2*z))*3.76*MWN2/Wfgasses;

% Caracteristics of the components 
% 0 means at 0 °C, a means at ambiant temperature, f means at t_exh

cpO20  = janaf('c','O2',273.15);
cpO2a  = janaf('c','O2',273.15 + Ta);
cpO2f  = janaf('c','O2',273.15 + t_exh);
cpN20  = janaf('c','N2',273.15);
cpN2a  = janaf('c','N2',273.15 + Ta);
cpN2f  = janaf('c','N2',273.15 + t_exh);
cpCO20 = janaf('c','CO2',273.15);
cpCO2f = janaf('c','CO2',273.15 + t_exh);
cpH2O0 = janaf('c','H2O',273.15);
cpH2Of = janaf('c','H2O',273.15 + t_exh);

% cp for ambiant air
cpAir_0a = 0.5 * (0.21*(cpO20 + cpO2a) + 0.79*(cpN20 + cpN2a));

% cp for flue gasses
cpFG = 0.5 * (CO2_frac*(cpCO20 + cpCO2f) + H2O_frac*(cpH2O0 + cpH2Of) ...
             + O2_frac*(cpO20  +  cpO2f) +  N2_frac*(cpN20  +  cpN2f));
         
epsilon = ((lambda*ma1 + 1)*cpFG*t_exh - lambda*ma1*cpAir_0a*Ta)/LHVCH4;

















end

