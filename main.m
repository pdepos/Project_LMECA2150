%% =======================================
% Project LMECA2150
%  by Philippe de Posson   5706-10-00
%     Jeremie Waltzing     xxxx-xx-xx
%% =======================================

Tmax = 540; % Max Steam Temperature
Pmax = 40;  % Max Steam Pressure
FH   = 6;   % Number of FeedHeaters 
RH   = 'off';   % Reheating on or off (2)

SimpleRankine   = RankineHirnOLD(540,40);
Steam_Cycle     = SteamCycle (Tmax, Pmax, FH, RH);

h     = SimpleRankine(:,4);
W_mcy = (h(3)-h(4)) - (h(2) - h(1));
Q1    = h(3)- h(2);

eta_cyclen = W_mcy / Q1;

disp('       States        t [°C]      p [bar]        x [-]    h [kJ/kg]   s [kJ/kgK]     e[kJ/kg] ')
%disp(SimpleRankine)
disp(Steam_Cycle)
