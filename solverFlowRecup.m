function [ma,mc,mg,lambda,T3,state] = solverFlowRecup(state,x,y,z,kmec,Pe,NTU,T5,fuel)

h1 = state{1}.h;
h2 = state{2}.h;
h4 = state{4}.h;
h5 = state{5}.h;

%%%% air characteristics %%%%
% 21% O2 et 79% N2
Ma = 0.21*32 + 0.79*28; % [g/mol]

mO2_a = 0.21*(32/Ma); %[kg_O2/kg_air]
mN2_a = 0.79*(28/Ma); %[kg_N2/kg_air]

%%%% Fuel characteristics %%%%

% Mfuel = x*12+y+z*16; % molar mass
if strcmp(fuel,'CH4')
    LHV = 50.1*10^3; %[kJ/kg]
elseif strcmp(fuel,'C12H23')
    LHV = 41.76*10^3; %[kJ/kg]
end
% DrH = 282.4*((2*z+y*x-y)/(2+y)) + ((2*x-2*z+y)/(2+y))*((1-y/4)*393.4 + (y/4)*802.4);
% LHV = (DrH/Mfuel)*10^3; %[kJ/kg] --> lower heating value

ma1 = ((32+3.76*28)*(x+(y-2*z)/4))/(12*x+y+16*z); % pouvoir comburivore [kg_air/kg_fuel]

%%%% Solver %%%%

T3 = (NTU/(1+NTU))*(T5-state{2}.T) + state{2}.T;

h0a = janaf('h','N2',273.15)*mN2_a + janaf('h','O2',273.15)*mO2_a;
state{3}.h = janaf('h','N2',T3)*mN2_a + janaf('h','O2',T3)*mO2_a - h0a;
h3 = state{3}.h;

ma = Pe/((1-kmec)*(1+(h4-h3)/abs(LHV-h4))*(h4-h5)-(1+kmec)*(h2-h1));
mc = ma*(h4-h3)/abs(LHV-h4);
mg = ma+mc;
lambda = ma/(ma1*mc);

end