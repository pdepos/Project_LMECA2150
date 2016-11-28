%% =======================================
% Project LMECA2150
%  by Philippe de Posson   5706-10-00
%     Jeremie Waltzing     xxxx-xx-xx
%% =======================================

Tmax = 565; % Max Steam Temperature
Tmin = 33;  % Temp at Condenser
Pmax = 170;  % Max Steam Pressure
FH   = 7;   % Number of FeedHeaters 
RH   = 'on';   % Reheating on or off

Steam_Cycle     = SteamCycle (Tmax, Tmin, Pmax, FH, RH);

disp('  States        t          p           x        h         s           e    ')
disp('              [�C]       [kPa]        [-]    [kJ/kg]   [kJ/kgK]    [kJ/kg] ')
disp('___________________________________________________________________________')
A = size(Steam_Cycle);
T = ones(1,A(1));
H = ones(1,A(1));
S = ones(1,A(1));
for i=1:A(1)
    fprintf('%8s %10.2f %10.1f %10.3f %10.1f %10.3f %10.1f\n',Steam_Cycle{i}.States,Steam_Cycle{i}.t,100*Steam_Cycle{i}.p, ...
        Steam_Cycle{i}.x,Steam_Cycle{i}.h,Steam_Cycle{i}.s,Steam_Cycle{i}.e);
   
    T(i) = Steam_Cycle{i}.t;
    H(i) = Steam_Cycle{i}.h;
    S(i) = Steam_Cycle{i}.s;
end

%Trying to display the interesting points in a graph. 

figure
subplot(2,1,1);
plot(S,T,'o');
subplot(2,1,2);
plot(S,H,'o');


