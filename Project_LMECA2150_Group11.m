%% =======================================
% Project LMECA2150
%  by Philippe de Posson   5706-10-00
%     Jeremie Waltzing     xxxx-xx-xx
%% =======================================

%% Input Data
%  Ici viendront les données internes et les input de l'interface graphique

ps3 = 310;   %[Bar]  Pression Steam point 3
ts3 = 565;   %[°C]   Temperature Steam point 3
ts5 = 565;
ts7 = 32;    %[°C]   Vérifier si unités OK pour calculs
ts6 = 32;
xs7 = 0;

eta_Si_HP = 0.92;
eta_Si_MP = 0.90;
eta_Si_BP = 0.88;

ns  = 23;    %       Number of states
% Fonction steam: Out=XSteam(‘function name’,In1,In2)

%% Steam Power Plant
t_s = ones(ns,1);
p_s = ones(ns,1);
x_s = ones(ns,1);
h_s = ones(ns,1);
s_s = ones(ns,1);
e_s = ones(ns,1);

p_s(3) = ps3;
t_s(3) = ts3;

%test: 
hs3 = XSteam('h_pt',ps3,ts3);
ss3 = XSteam('s_pt',ps3,ts3);

ps7 = XSteam('psat_t',ts7);
hs7 = XSteam('hL_p',ps7);
ss7 = XSteam('sL_p',ps7);
