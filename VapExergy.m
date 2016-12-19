function [stateV] = VapExergy(stateV)

n = length(stateV);

T0 = 15; % [°C]
p0 = 1; % [bar]
h0 = XSteam('h_pT',p0,T0);
s0 = XSteam('s_pT',p0,T0);

for i = 1:n
   stateV{i}.e = (stateV{i}.h-h0) - (T0+273.15)*(stateV{i}.s-s0); 
end

end