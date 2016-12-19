%% =======================================
% Project LMECA2150
%  by Philippe de Posson   5706-10-00
%     Jeremie Waltzing     xxxx-xx-xx
%% =======================================
close all

Ta   = 15;         % Ambient Air Temperature
Triv = 15;         % River Temperature   
Texh = 120;        % Temperature exhaust flue gasses
Tmax = 560;        % Max Steam Temperature
Tmin = Triv + 18;  % Temp at Condenser
Pmax = 50;         % Max Steam Pressure
P_el = 35e3;       % Electric Power of the Plant in kW
FH   = 0;          % Number of FeedHeaters 
RH   = 'off';      % Reheating on or off

% Fuel Type: CxHyOz
x = 1;
y = 4;
z = 0;
lambda = 1.05;

Steam_Cycle = SCStates (Tmax, Tmin, Pmax, FH, RH);
[ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] = SCAnalysis (Steam_Cycle, P_el, FH, RH, x, y ,z, Texh, Ta, lambda);

if FH == 0
    A = length(Steam_Cycle);
else 
    A = length(Steam_Cycle) - 1;
end

skipped = 0;
States = cell(A,7);
for i = 1:A(1)
    if Steam_Cycle{i}.p == 0
        skipped = 1;
    end
    if skipped == 0 
        States{i,1} = Steam_Cycle{i}.States;
        States{i,2} = Steam_Cycle{i}.t;
        States{i,3} = 100*Steam_Cycle{i}.p;
        States{i,4} = Steam_Cycle{i}.x;
        States{i,5} = Steam_Cycle{i}.h;
        States{i,6} = Steam_Cycle{i}.s;
        States{i,7} = Steam_Cycle{i}.e;  
    else
        States{i,1} = Steam_Cycle{i+1}.States;
        States{i,2} = Steam_Cycle{i+1}.t;
        States{i,3} = 100*Steam_Cycle{i+1}.p;
        States{i,4} = Steam_Cycle{i+1}.x;
        States{i,5} = Steam_Cycle{i+1}.h;
        States{i,6} = Steam_Cycle{i+1}.s;
        States{i,7} = Steam_Cycle{i+1}.e;      
    end
end

% 
% disp('  States        t          p           x        h         s           e    ')
% disp('              [°C]       [kPa]        [-]    [kJ/kg]   [kJ/kgK]    [kJ/kg] ')
% disp('___________________________________________________________________________')
% for i=1:A(1)
%     fprintf('%8s %10.2f %10.1f %10.3f %10.1f %10.3f %10.1f\n',Steam_Cycle{i}.States,Steam_Cycle{i}.t,100*Steam_Cycle{i}.p, ...
%         Steam_Cycle{i}.x,Steam_Cycle{i}.h,Steam_Cycle{i}.s,Steam_Cycle{i}.e);
% end

% str = cell(length(EnergyDistribution)-1,1)
% for i = 1:length(EnergyDistribution)-1
% str{i} = fprintf(formatSpec,str1,value)

formatSpec = '%s %.2f';
str1   = EnergyDistribution{1,1};
value1 = EnergyDistribution{2,1}/1000;
label1 = sprintf(formatSpec,str1,value1);
str2   = EnergyDistribution{1,2};
value2 = EnergyDistribution{2,2}/1000;
label2 = sprintf(formatSpec,str2,value2);
str3   = EnergyDistribution{1,3};
value3 = EnergyDistribution{2,3}/1000;
label3 = sprintf(formatSpec,str3,value3);
str4   = EnergyDistribution{1,4};
value4 = EnergyDistribution{2,4}/1000;
label4 = sprintf(formatSpec,str4,value4);
   


Labels = {label1 label2 label3 label4};
figure 
%subplot(3,1,2)
h = pie( [EnergyDistribution{2,1:4}],Labels);%EnergyDistribution(2,1:4));%Labels);

hText = findobj(h,'Type','text'); % text object handles
percentValues = get(hText,'String'); % percent values

oldExtents_cell = get(hText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell); % numeric array


newExtents_cell = get(hText,'Extent'); % cell array
newExtents = cell2mat(newExtents_cell); % numeric array
width_change = newExtents(:,4)-oldExtents(:,4);

signValues = sign(oldExtents(:,1));
offset = 0.9*signValues;%.*(width_change/2);

textPositions_cell = get(hText,{'Position'}); % cell array
textPositions = cell2mat(textPositions_cell); % numeric array
textPositions(:,1) = textPositions(:,1) + offset; % add offset

hText(1).Position = textPositions(1,:);
hText(2).Position = textPositions(2,:);
hText(3).Position = textPositions(3,:);
hText(4).Position = textPositions(4,:);



% embiggenby = 0; % <-enter a value here; it's a percent. 
% % Make a plot: 
% %subplot(2,1,2) 
% pie(rand(10,1))
% % Get position of the plot: 
% pos = get(gca,'outerposition'); 
% % Change axis position: 
% pos(2) = pos(2) - embiggenby/2; 
% pos(3) = pos(3)+ embiggenby 
% % Set new axis position: 
% set(gca,'outerposition',pos)








%SCPlot(Steam_Cycle,RH,FH);
% figure
% subplot(2,1,1);
% pie( [EnergyDistribution{2,1:4}],EnergyDistribution(1,1:4)');
% subplot(2,1,2);
% pie( [ExergyDistribution{2,:}],ExergyDistribution(1,:)');
% % 

% % 
%   FlowRates
%  eta_en
%  EnergyDistribution
%   eta_ex
%   ExergyDistribution

% PrimaryEn = P_el/eta_en{2,4};
% PrimaryEnSum = 0;
% for i = 1:length(EnergyDistribution) - 1
%     PrimaryEnSum = PrimaryEnSum + EnergyDistribution{2,i};
% end
% EnError = (PrimaryEn - PrimaryEnSum);
% PrimaryEx = P_el/eta_ex{2,2};
% PrimaryExSum = 0;
% for i = 1:length(ExergyDistribution)
%     PrimaryExSum = PrimaryExSum + ExergyDistribution{2,i};
% end
% ExError = (PrimaryEx - PrimaryExSum);
% 
% fprintf('Primary Energy Flux is %.2f MW, PrimaryEnSum is %.2f MW, EnError is %f MW\n',PrimaryEn/1000, PrimaryEnSum/1000, EnError/1000);
% fprintf('Primary Exergy Flux is %.2f MW, PrimaryExSum is %.2f MW, ExError is %f MW\n',PrimaryEx/1000, PrimaryExSum/1000, ExError/1000);