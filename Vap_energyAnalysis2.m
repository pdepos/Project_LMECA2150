function [Pcond,PfmecV,etaMecV,etaCyclenV,etaTotenV] = Vap_energyAnalysis2(stateV,mvHP,mvLP)

kmec = 0.015;

%%%% losses at the condenser %%%%
Pcond = (mvHP+mvLP)*(stateV{10}.h - stateV{1}.h); %[kW]

%%%% mechanical losses %%%%
PmT = mvHP*(stateV{8}.h-stateV{9}.h) + (mvHP+mvLP)*(stateV{9}.h-stateV{10}.h); % [kW]
PfmecV = kmec*PmT; % [kW]

%%%% Heat gained in the exchanger %%%%
Q = mvHP*(stateV{8}.h-stateV{5}.h) + mvLP*(stateV{9}.h-stateV{3}.h) +...
    (mvHP+mvLP)*(stateV{3}.h-stateV{2}.h); % [kW]

%%%% Efficiencies %%%%
etaMecV = 1 - kmec;
etaCyclenV = PmT/Q;
etaTotenV = etaMecV*etaCyclenV;
end