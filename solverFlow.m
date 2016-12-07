function [ma,mc,mg,lambda] = solverFlow(state,x,y,z,kmec,Pe,fuel)


h1 = state{1}.h;
h2 = state{2}.h;
h3 = state{3}.h;
h4 = state{4}.h;

%%%% Fuel characteristics %%%%

Mfuel = x*12+y+z*16; % molar mass
if strcmp(fuel,'CH4')
    LHV = 50.1*10^3; %[kJ/kg]
elseif strcmp(fuel,'C12H23')
    LHV = 41.76*10^3; %[kJ/kg]
end
% DrH = 282.4*((2*z+y*x-y)/(2+y)) + ((2*x-2*z+y)/(2+y))*((1-y/4)*393.4 + (y/4)*802.4);
% LHV = (DrH/Mfuel)*10^3; %[kJ/kg] --> lower heating value

ma1 = ((32+3.76*28)*(x+(y-2*z)/4))/(12*x+y+16*z); % pouvoir comburivore [kg_air/kg_fuel]

%%%% Solver %%%%

ma = Pe/((1-kmec)*(1+(h3-h2)/abs(LHV-h3))*(h3-h4)-(1+kmec)*(h2-h1));
mc = ma*(h3-h2)/abs(LHV-h3);
mg = ma+mc;
lambda = ma/(ma1*mc);
end