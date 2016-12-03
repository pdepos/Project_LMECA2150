%% =======================================
% Project LMECA2150
%  by Philippe de Posson   5706-10-00
%     Jeremie Waltzing     xxxx-xx-xx
%% =======================================

Tmax  = 520; % Max Steam Temperature
Tmin  = 33;  % Temp at Condenser
T_exh = 120; % Temperature exhaust flue gasses
Ta    = 15;
Pmax  = 40; % Max Steam Pressure
P_el_LP = 35e3;  % Electric Power of the LPT in kW
P_el_HP = 100e3;  % Electric Power of the HPT in kW
FH   = 0;   % Number of FeedHeaters 
RH   = 'off';% Reheating on or off

%Fuel Type: CxHyOz
x = 1;
y = 4;
z = 0;
lambda = 1.05;

Steam_Cycle = SteamCycle (Tmax, Tmin, Pmax, FH, RH);
[ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution, T_exh ] = EnergyAnalysisSC (Steam_Cycle, P_el_LP, P_el_HP, FH, RH, x, y ,z, T_exh, Ta, lambda);


disp('  States        t          p           x        h         s           e    ')
disp('              [°C]       [kPa]        [-]    [kJ/kg]   [kJ/kgK]    [kJ/kg] ')
disp('___________________________________________________________________________')
A = size(Steam_Cycle);
T_Cycle = ones(1,A(1));
H_Cycle = ones(1,A(1));
S_Cycle = ones(1,A(1));
for i=1:A(1)
    fprintf('%8s %10.2f %10.1f %10.3f %10.1f %10.3f %10.1f\n',Steam_Cycle{i}.States,Steam_Cycle{i}.t,100*Steam_Cycle{i}.p, ...
        Steam_Cycle{i}.x,Steam_Cycle{i}.h,Steam_Cycle{i}.s,Steam_Cycle{i}.e);
   
    T_Cycle(i) = Steam_Cycle{i}.t;
    H_Cycle(i) = Steam_Cycle{i}.h;
    S_Cycle(i) = Steam_Cycle{i}.s;
end

T  = linspace(0.1,373.945813,2000);
HL = zeros(1,length(T));
SL = zeros(1,length(T));
HV = zeros(1,length(T));
SV = zeros(1,length(T));

for i = 1:length(T)
    HL(i) = XSteam('hL_T',T(i));
    SL(i) = XSteam('sL_T',T(i));
    HV(i) = XSteam('hV_T',T(i));
    SV(i) = XSteam('sV_T',T(i));
end

figure
subplot(2,1,1);
pie( [EnergyDistribution{2,:}],EnergyDistribution(1,:)');
subplot(2,1,2);
pie( [ExergyDistribution{2,:}],ExergyDistribution(1,:)');

fprintf('Temperature exhaust gasses: %.0f °C.',T_exh);
FlowRates
eta_en
EnergyDistribution
eta_ex
ExergyDistribution

figure
subplot(2,1,1);
plot(S_Cycle,T_Cycle,'or',SL,T,'-b',SV,T,'-b');
subplot(2,1,2);
plot(S_Cycle,H_Cycle,'or',SL,HL,'-b',SV,HV,'-b');


