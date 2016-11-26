function [state] = GasTurbineExergy(state)
% function computing the exergy of each state of a Gas Turbine


n = length(state);

Ma = 0.21*32 + 0.79*28; % [g/mol]
mO2_a = 0.21*(32/Ma); %[kg_O2/kg_air]
mN2_a = 0.79*(28/Ma); %[kg_N2/kg_air]

T0 = 273.15 + 15; %[K]

h0a = janaf('h','N2',273.15)*mN2_a + janaf('h','O2',273.15)*mO2_a;
s0a = janaf('s','N2',273.15)*mN2_a + janaf('s','O2',273.15)*mO2_a;

h0 = janaf('h','N2',273.15 + 15)*mN2_a + janaf('h','O2',273.15 + 15)*mO2_a - h0a; %[KJ/kg]
s0 = janaf('s','N2',273.15 + 15)*mN2_a + janaf('s','O2',273.15 + 15)*mO2_a - s0a; %[KJ/(kg*K)]

for i=1:n
    state{i}.e = (state{i}.h-h0)-T0*(state{i}.s-s0);
end

end