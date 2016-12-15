function ec = fuel_exergy(x,y,z)
% function calculating the exergy of the fuel CxHyOz in [kJ/kg]

%%%% Data %%%%
x = 12;
y = 23;
z = 0;
% molar mass [kg/kmol]
M_H = 1;
M_C = 12;
M_O = 16;

M_H2 = 2*M_H;
M_CO2 = M_C + 2*M_O;
M_O2 = 2*M_O;
M_H2O = M_H + 2*M_O;
M_CO = M_C + M_O;
M_CH4 = M_C + 4*M_H;

% average Cp between 273 K and 288 K [kJ/kg*K]
Cp_H2 = 28.7/M_H2;
Cp_O2 = 29.3/M_O2;
Cp_CO2 = 36.5/M_CO2;
Cp_H2O_liq = 75.5/M_H2O;
Cp_C = 10.4/M_C;
Cp_CO = 29.1/M_CO;
Cp_CH4 = 14.4/M_CH4;

% entropy at reference state [kJ/kg]
S_O2 = 202.8/M_O2;
S_CO2 = 210.4/M_CO2;
S_H2O_liq = 69.5/M_H2O;
S_C = 4.78/M_C;
S_CO = 194.8/M_CO;
S_CH4 = 183.1/M_CH4;

% gas constant
R = 8.314472;
R_O2 = R/M_O2;
R_CO2 = R/M_CO2;

%%%% Fuel characteristics %%%%
Mfuel = x*12+y+z*16; % molar mass
DrH = 282.4*((2*z+y*x-y)/(2+y)) + ((2*x-2*z+y)/(2+y))*((1-y/4)*393.4 + (y/4)*802.4);
LHV = (DrH/Mfuel)*10^3 %[kJ/kg] --> lower heating value

Hv_H2O = 2511; % heat of vaporization [kJ/kg]
HHV = LHV + Hv_H2O*(y/2)*(18/Mfuel); % higher heating value in [kJ/kg]

Cp_fuel = ((2*z+y*(x-1))/(2+y))*Cp_CO*(M_CO/Mfuel) + ((y*(1+z-x))/(2+y))*Cp_H2O_liq*(M_H2O/Mfuel) + ...
    ((2*x-2*z+y)/(2+y))*((1-y/4)*Cp_C*(M_C/Mfuel) + y/4*Cp_CH4*(M_CH4/Mfuel)); % [kJ/(kg*K)]

S_fuel = ((2*z+y*(x-1))/(2+y))*S_CO*(M_CO/Mfuel) + ((y*(1+z-x))/(2+y))*S_H2O_liq*(M_H2O/Mfuel) + ...
    ((2*x-2*z+y)/(2+y))*((1-y/4)*S_C*(M_C/Mfuel) + y/4*S_CH4*(M_CH4/Mfuel));% [kJ/(kg*K)]

%%%% calculation of mass fractions %%%%
m_O2 = (x+(y-2*z)/4)*(1/Mfuel)*M_O2;
m_CO2 = x*(1/Mfuel)*M_CO2;
m_H2O = (y/2)*(1/Mfuel)*M_H2O;

%%%% calculation of exergy %%%%
ec = HHV + 15*(Cp_fuel + m_O2*Cp_O2 - m_CO2*Cp_CO2 - m_H2O*Cp_H2O_liq)...
    - 288*(S_fuel + Cp_fuel*log(288/273))...
    - 288*m_O2*(S_O2 + Cp_O2*log(288/273) - R_O2*log(20.63/100))...
    + 288*m_CO2*(S_CO2 + Cp_CO2*log(288/273) - R_CO2*log(0.03/100))...
    + 288*m_H2O*(S_H2O_liq + Cp_H2O_liq*log(288/273))

end