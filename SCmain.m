%% =======================================
% Project LMECA2150
%  by Philippe de Posson   5706-10-00
%     Jeremie Waltzing     xxxx-xx-xx
%% =======================================
close all

Texh  = 120;        % Temperature exhaust flue gasses
Ta    = 15;         % Ambient Air Temperature
Tmax  = 520;        % Max Steam Temperature
Triv  = 15;         % River Temperature    
Tmin  = Triv + 18;  % Temp at Condenser
Pmax  = 40;         % Max Steam Pressure
P_el_LP = 35e3;     % Electric Power of the LPT in kW
P_el_HP = 65e3;     % Electric Power of the HPT in kW
FH   = 7;           % Number of FeedHeaters 
RH   = 'off';       % Reheating on or off

% Fuel Type: CxHyOz
x = 1;
y = 4;
z = 0;
lambda = 1.05;

Steam_Cycle = SCStates (Tmax, Tmin, Pmax, FH, RH);
[ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] = SCAnalysis (Steam_Cycle, P_el_LP, P_el_HP, FH, RH, x, y ,z, Texh, Ta, lambda);


% disp('  States        t          p           x        h         s           e    ')
% disp('              [°C]       [kPa]        [-]    [kJ/kg]   [kJ/kgK]    [kJ/kg] ')
% disp('___________________________________________________________________________')
% A = size(Steam_Cycle);
% for i=1:A(1)
%     fprintf('%8s %10.2f %10.1f %10.3f %10.1f %10.3f %10.1f\n',Steam_Cycle{i}.States,Steam_Cycle{i}.t,100*Steam_Cycle{i}.p, ...
%         Steam_Cycle{i}.x,Steam_Cycle{i}.h,Steam_Cycle{i}.s,Steam_Cycle{i}.e);
% end


%SCPlot(Steam_Cycle,RH,FH);



% figure
% subplot(2,1,1);
% pie( [EnergyDistribution{2,:}],EnergyDistribution(1,:)');
% subplot(2,1,2);
% pie( [ExergyDistribution{2,:}],ExergyDistribution(1,:)');
% 
FlowRates
% eta_en
% EnergyDistribution
% eta_ex
ExergyDistribution

PrimaryExtest = 0;
for i = 1:length(ExergyDistribution)
    PrimaryExtest = PrimaryExtest + ExergyDistribution{2,i};
end

PrimaryEx = P_el_LP/eta_ex{2,2};
fprintf('Primary Exergy Flux is %.2f MW, Primary Extest is %.2f MW \n',PrimaryEx/1000, PrimaryExtest/1000);


% PrimaryEnergy = EnergyDistribution{2,1}/eta_en{2,length(eta_en)}
% PrimaryExergy = ExergyDistribution{2,1}/eta_ex{2,2}




