%% =======================================
% Project LMECA2150
%  by Philippe de Posson   5706-10-00
%     Jeremie Waltzing     xxxx-xx-xx
%% =======================================

SimpleRankine = RankineHirn(40,540);

h     = SimpleRankine(:,4);
W_mcy = (h(3)-h(4)) - (h(2) - h(1));
Q1    = h(3)- h(2);

eta_cyclen = W_mcy / Q1;

disp('       t [°C]      p [bar]        x [-]    h [kJ/kg]   s [kJ/kgK]     e[kJ/kg] ')
disp(SimpleRankine)