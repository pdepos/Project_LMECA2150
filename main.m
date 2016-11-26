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

disp('  States        t          p           x        h         s           e    ')
disp('              [°C]       [kPa]        [-]    [kJ/kg]   [kJ/kgK]    [kJ/kg] ')
disp('___________________________________________________________________________')
A = size(Steam_Cycle);
T = ones(1,A(1));
S = ones(1,A(1));
for i=1:A(1)
    fprintf('%8s %10.2f %10.1f %10.3f %10.1f %10.3f %10.1f\n',Steam_Cycle{i}.States,Steam_Cycle{i}.t,100*Steam_Cycle{i}.p, ...
        Steam_Cycle{i}.x,Steam_Cycle{i}.h,Steam_Cycle{i}.s,Steam_Cycle{i}.e);
   
    T(i) = Steam_Cycle{i}.t;
    S(i) = Steam_Cycle{i}.s;
end

%Trying to display the interesting points in a graph. 

plot(S,T,'o');