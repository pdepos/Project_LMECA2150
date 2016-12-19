function [Pcond,PfmecV,etaMecV,etaCyclenV,etaTotenV] = Vap_energyAnalysis3(stateV,mvHP,mvMP,mvLP)

kmec = 0.015;

%%%% losses at the condenser %%%%
Pcond = (mvHP+mvMP+mvLP)*(stateV{20}.h - stateV{1}.h); %[kW]

%%%% mechanical losses %%%%
PmT = mvHP*(stateV{14}.h-stateV{15}.h) + (mvHP+mvMP)*(stateV{17}.h-stateV{18}.h)...
    + (mvHP+mvMP+mvLP)*(stateV{19}.h-stateV{20}.h); % [kW]
PfmecV = kmec*PmT; % [kW]

%%%% Heat gained in the exchanger %%%%
Q = mvHP*(stateV{14}.h-stateV{9}.h) + mvMP*(stateV{13}.h-stateV{8}.h) +...
    mvLP*(stateV{6}.h-stateV{3}.h) + (mvHP+mvMP+mvLP)*(stateV{3}.h-stateV{2}.h) + ...
    (mvHP+mvMP)*(stateV{17}.h-stateV{16}.h+stateV{8}.h-stateV{7}.h); % [kW]

%%%% Efficiencies %%%%
etaMecV = 1 - kmec;
etaCyclenV = PmT/Q;
etaTotenV = etaMecV*etaCyclenV;

end