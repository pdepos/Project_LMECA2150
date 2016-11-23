function [state] = GasTurbineExergy(state)
% function computing the exergy of each state of a Gas Turbine


n = length(state);

T0 = 273.15 + 15; %[K]
h0 = 15.1; %[KJ/kg]
s0 = 0.054; %[KJ/(kg*K)]

for i=1:n
    state{i}.e = (state{i}.h-h0)-T0*(state{i}.s-s0);
end

end