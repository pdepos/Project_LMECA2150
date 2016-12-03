function [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution, T_exh] = EnergyAnalysisSC( state, P_el_LP, P_el_HP, FH, RH, x, y, z, t_exh, Ta, lambda )
%Steam Cycle
%   Input Arguments: 
%   - Max Steam Pressure [bar]
%   - Max Temperature [°C]
%   - Number of Feed Heaters [min value is 0]
%   - Reheating on or off
%   =======================================================================

%%
% Global Variables
MWO2   = 31.998;
MWN2   = 28.014;
MWCO2  = 44.008;
MWH2O  = 18.01494;
MWH2   = 2.01594;
MWair  = 0.21*MWO2 + 0.79*MWN2;
MWC    = MWCO2 - MWO2;
MWCH4  = MWC + 2*MWH2;
MWCH18 = MWC + 1.8 / 2 * MWH2;

if (x == 1) && (y == 4) && (z == 0)
    LHV  = 50150; %kJ/kg
    HHV  = 32780; %kJ/kg
    ec   = 52215; %kJ/kg
    cp15 = 35.3 / MWCH4;
    f    = ec / LHV; % Exergy of methane / LHV methane
elseif (x == 12) && (y == 23) && (z == 0)
    LHV  = 42900; %LHV of CH1.8 which will approximate Diesel, considering that the LHV of heavy fuels varies less as that for gasses. 
    HHV  = 45790;
    ec   = 45710;
    cp15 = 14.4 / MWCH18;
    f    = ec / LHV;
else
    warning('Unrecognized type of fuel. Default value: CH4')
    LHV  = 50150; %kJ/kg
    HHV  = 32780; %kJ/kg
    ec   = 52215; %kJ/kg
    cp15 = 35.3 / MWCH4;
    f    = ec / LHV; % Exergy of methane / LHV methane
end
if lambda < 1
    warning('Lambda cannot be under 1 in our study. Default value: lambda = 1')
    lambda = 1;
end

eta_mec = 0.98;
ma1     = ((MWO2+3.76*MWN2)*(x + (y - 2*z)/4))/(12*x + y + 16*z);

switch RH
    case 'off' 
        r = 0;
        if FH == 0
            [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution, T_exh ] ...
                = RankineHirnEnergy(state, P_el_LP, t_exh);
        elseif FH > 0
            [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] ...
                = BasicFeedHeatingEnergy(state, P_el_LP, FH, t_exh, r);
        else
            warning('Negative number of feedheaters not allowed')
            EnergyDistribution = 0;
            FlowRates = 0;
        end
    case 'on'
        r = 1;
        if FH == 0
            [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] ...
                = BasicReHeating(T_max, T_min, P_max);
        elseif FH > 0
            [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] ...
                = ReHeatingAndFHEnergy(state, P_el_LP, P_el_HP, FH, r);
        else
            warning('Negative number of feedheaters not allowed')
            EnergyDistribution = 0;
            FlowRates = 0;
        end
    otherwise
        warning('Unexpected reheating entry.')
        EnergyDistribution = 0;
        FlowRates = 0;
end


%%
function [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution, T_exh] = RankineHirnEnergy(state, P_el, t_exh)

ind1   = 1;
ind2   = 2;
ind3   = 5;
ind4   = 6;

% Vapor Flow Rate
m_dot_v = P_el / ((state{ind3}.h - state{ind4}.h)*eta_mec);

% Energetic Efficiencies
eta_cyclen = ((state{ind3}.h - state{ind4}.h) - abs(state{ind2}.h - state{ind1}.h)) / ...
             (state{ind3}.h - state{ind2}.h);
eta_gen = Efficiency_Gen(t_exh);
eta_toten = eta_mec*eta_gen*eta_cyclen;


% Energetic Losses
Q_gen = m_dot_v * (state{ind3}.h - state{ind2}.h); % Power recieved by the fluid at generator
Q_comb = Q_gen / eta_gen;                          % Primary Power needed from the fuel
ConEnLosses  = m_dot_v * (state{ind4}.h - state{ind1}.h);
MecEnLosses  = (1 - eta_mec) * P_el;
GenEnLosses  = (1 - eta_gen) * Q_comb;

% Combustible Flow Rate
m_dot_c = Q_comb / LHV;

% Flue Gasses Flow Rate
m_dot_f = m_dot_c * (1 + lambda*ma1);

% Exergetic Efficiencies
cpair = cpair(Ta);   % Verified OK
ha    =  cpair*Ta;   % Verified OK
% hatest= EnthalpyAir(Ta) - EnthalpyAir(0); % Gives about the same as cpair*ta

hc    =   15*cp15;   % Verified OK
hfg_comb    = (LHV + hc + lambda*ma1*ha)/(lambda*ma1 + 1);                     % Verified OK
hfg_exitboiler = hfg_comb - m_dot_v/m_dot_f * (state{ind3}.h - state{ind2}.h); % Verified OK

Tcomb        = TemperatureFG (hfg_comb);             % is ok
T_exitboiler = TemperatureFG (hfg_exitboiler);       % is ok

efg_comb       = ExergyFlueGasses(hfg_comb, Tcomb);                % is ok
efg_exitboiler = ExergyFlueGasses(hfg_exitboiler, T_exitboiler);   % is ok
er             = 0;                                  % Exergy Reagents = 0 at T_reagents = 15°C

eta_totex   = eta_toten / f;
eta_combex  = m_dot_f*(efg_comb - er) / (m_dot_c * ec);
eta_chimnex = m_dot_f*(efg_comb - efg_exitboiler)/(m_dot_f*(efg_comb - er));
eta_transex = m_dot_v*(state{ind3}.e - state{ind2}.e)/(m_dot_f*(efg_comb - efg_exitboiler));
eta_gex     = eta_combex*eta_chimnex*eta_transex;
%eta_gex    = m_dot_v*(state{ind3}.e - state{ind2}.e) / (m_dot_c * ec);
eta_rotex   = ((state{ind3}.h - state{ind4}.h) - (state{ind2}.h - state{ind1}.h)) /...
              ((state{ind3}.e - state{ind4}.e) - (state{ind2}.e - state{ind1}.e));
eta_cyclex  = ((state{ind3}.h - state{ind4}.h) - (state{ind2}.h - state{ind1}.h)) /...
              (state{ind3}.e - state{ind2}.e);
          
% Exergetic Losses
MecExLosses = MecEnLosses;
ConExLosses = m_dot_v * (state{ind4}.e - state{ind1}.e);
RotExLosses = (1 -   eta_rotex) * m_dot_v * ((state{ind3}.e - state{ind4}.e) - (state{ind2}.e - state{ind1}.e));
TraExLosses = (1 - eta_transex) * m_dot_f * (efg_comb - efg_exitboiler);
StaExLosses = (1 - eta_chimnex) * m_dot_f * (efg_comb - er);
ComExLosses = (1 -  eta_combex) * m_dot_c * ec;
          
% Outputs
FlowRates = {'m°v' 'm°c' 'm°f'; m_dot_v m_dot_c m_dot_f};
eta_en    = {'eta_mec', 'eta_gen', 'eta_cyclen', 'eta_toten';...
            eta_mec, eta_gen, eta_cyclen, eta_toten};
EnergyDistribution = ...
    {'Effective Power', 'Mechanical Losses', 'Condensor Losses', 'Generator Losses';...
    P_el, MecEnLosses, ConEnLosses, GenEnLosses};
eta_ex = ...
    {'eta_mec' 'eta_totex' 'eta_gex' 'eta_combex' 'eta_chimnex' 'eta_transex' 'eta_rotex' 'eta_cyclex';...
    eta_mec eta_totex eta_gex eta_combex eta_chimnex eta_transex eta_rotex eta_cyclex};
ExergyDistribution = ...
    {'Effective Power' 'Mechanical Losses' 'Condensor Losses' 'Irreversibilities in the Turbine and Pumps'...
    'Heat Transfer Irreversibility in the Steam Generator' 'Stack Losses' 'Combustion Irreversibility';...
    P_el MecExLosses ConExLosses RotExLosses TraExLosses StaExLosses ComExLosses};

T_exh = T_exitboiler;
end
%%
function [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] = BasicFeedHeatingEnergy(state, P_el, fh, x, y, z, t_exh, Ta, Lambda, r)

    
alpha = 10;
beta  = 4;  
gamma = 4;

ind1   = 1;
ind2   = 2;
ind2i  = 3;
ind2ii = 4;
ind3   = 5;
ind4   = 6;
ind5   = 7;
ind6   = 8;
ind7   = 9;
ind8   = 10;
ind41  = 11;

n = alpha + beta * fh;

% Left Hand Side of the flow rates system
A = zeros(fh,1);
% Right Hand Side of the flow rates system
B = zeros(fh);

position_bache = 1;               %References the roman number of the schematic
if fh == 1
    A(1) = state{ind1}.h - state{ind6}.h;
    B(1) = state{ind41}.h - state{ind7}.h - A(1);
elseif fh > gamma                     %Finding the bache if the number of bleeds is higher than gamma
    while state{index_bleed_bache}.p < 4.6
        position_bache = position_bache + 1;
    end
    for i = 1:fh
        ind6i     = alpha + beta*i - 2;
        ind6iplus = alpha + beta*(i + 1) - 2;
        ind4i     = alpha + beta*(i-1) + 1; 
        ind7i     = alpha + beta*i - 1;
        ind8iplus = alpha + beta*(i + 1);

        if i == 1
            A(i) = state{ind6iplus}.h - state{ind6}.h;
        elseif i == position_bache
            A(i) = state{ind7i}.h - state{ind6i}.h;
        elseif i == fh
            A(i) = state{ind1}.h - state{ind6i}.h;
        else A(i) = state{ind6iplus}.h - state{ind6i}.h;
        end

        for j = 1:fh
            % First Line
            if (i == 1) && (j == 1)
                B(i,j) = (state{ind41}.h - state{ind7}.h) - A(i);
            elseif (i == 1) && (j < position_bache) 
                B(i,j) = (state{ind8iplus}.h - state{ind7}.h) - A(i);

            % Coeffs above diagonal before the bache
            elseif (i > 1) && (j < position_bache) && (j > i)
                B(i,j) = (state{ind8iplus}.h - state{ind7i}.h) - A(i);

            % The flow of the bleeds occuring behind the bache are not part of
            % the feedwater before the bache
            elseif (i < position_bache) && (j >= position_bache)
                B(i,j) = 0;

            % Coeffs under the array diagonal
            elseif (j < i)
                B(i,j) = - A(i);

            % Coeffs on the diagonal except at the bache
            elseif (i > 1) && (i == j) && (i ~= position_bache)
                B(i,j) = (state{ind4i}.h - state{ind7i}.h) - A(i);

            % Coeff on the diagonal at the position of the bache and behind
            elseif (i == position_bache) && (j == position_bache)
                B(i,j) = state{ind4i}.h - state{ind7i}.h;       
            elseif (i == position_bache) && (j > position_bache)
                B(i,j) = state{ind8iplus}.h - state{ind7i}.h;
            elseif (i > position_bache) && ( j > i)
                B(i,j) = (state{ind8iplus}.h - state{ind7i}.h) - A(i);
            else
                B(i,j) = NaN;
            end
        end
    end
else
    for i = 1:fh
        ind6i     = alpha + beta*i - 2;
        ind6iplus = alpha + beta*(i + 1) - 2;
        ind4i     = alpha + beta*(i-1) + 1; 
        ind7i     = alpha + beta*i - 1;
        ind8iplus = alpha + beta*(i + 1);

        if i == 1
            A(i) = state{ind6iplus}.h - state{ind6}.h;
        elseif i == fh
            A(i) = state{ind1}.h - state{ind6i}.h;
        else A(i) = state{ind6iplus}.h - state{ind6i}.h;
        end

        for j = 1:fh
            % First Line
            if (i == 1) && (j == 1)
                B(i,j) = (state{ind41}.h - state{ind7}.h) - A(i);
            elseif (i == 1) && (j > 1)
                B(i,j) = (state{ind8iplus}.h - state{ind7}.h) - A(i);

            % Coeffs above the diagonal
            elseif (i > 1) && (j > i)
                B(i,j) = (state{ind8iplus}.h - state{ind7i}.h) - A(i);

            % Coeffs under the diagonal
            elseif (j < i)
                B(i,j) = - A(i);

            % Coeffs on the diagonal 
            elseif (i > 1) && (i == j)
                B(i,j) = (state{ind4i}.h - state{ind7i}.h) - A(i);
                
            % To detect errors
            else
                B(i,j) = NaN;
            end
        end
    end
end

X = B\A
w_turbine = state{ind3}.h - state{ind4}.h;
for i = 1:fh
    ind4i     = alpha + beta*(i-1) + 1; 
    w_turbine = w_turbine + X(i)*(state{ind3}.h - state{ind4i}.h);  
end

P_meca = P_el/eta_mec;
main_flowrate = P_meca / w_turbine;

FR = ones(1,fh+1);
FR(1) = main_flowrate;
for i = 2:fh + 1
    FR(i) = FR(1)*X(i - 1);
end


FlowRates = FR
end
%%
function [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] = ReHeatingAndFHEnergy(state, P_el_LP, P_el_HP, fh, r);
alpha = 12;                %Number of fixed points of the schematic
beta  = 4;                 %Number of points per state: eg: beta = 3 means 4.1 6.1 7.1 or 4.2 6.2 7.2 
                   %It's in case we need to add one to take
                   %something more into account. eg before and
                   %after isenthalpic valves
alpha_FH_base = 10; 

ind1    = 1;
ind2    = 2;
ind2i   = 3;
ind2ii  = 4;
ind3HP  = 5;
ind4HP  = 6;
ind3LP  = 7;
ind4LP  = 8;
ind5    = 9;
ind6    = 10;
ind7    = 11;
ind8    = 12;

n = alpha + beta * fh + 4;  

P_meca_HP = P_el_HP / eta_mec;
W_HPT = state{ind3HP}.h - state{ind4HP}.h;
main_flowrate = P_meca_HP / W_HPT;


FlowRates = main_flowrate
end
%%
function [ eta_gen ] = Efficiency_Gen(t_exh)
%UNTITLED Summary of this function goes here
%   Input Args:
%       - CxHyOz: Fuel type
%       - t_exh: Temperature of the exhaust gasses
%       - lambda: Excess Air Coefficient, cannot be smaller than 1

cpAir_Ta = cpair(Ta);
cpFG = cpfg (t_exh);

epsilon = ((lambda*ma1 + 1)*cpFG*t_exh - lambda*ma1*cpAir_Ta*Ta)/LHV;

eta_gen = 1 - epsilon - 0.01;
end

%%
function [ exergyFG ] = ExergyFlueGasses(hf, tfg)
% Function Calculating the Exergy of the Flue Gasses at the combustion
% temperature. 

% % Finding the temperature of combustion
% tfg0 = 1500;                        % First guess for T combustion
% cpfg0 = cpfg (tfg0);   % First guess for cp flue gasses
% 
% tfgprev = tfg0;
% tfgnew  = hf / cpfg0;
% while abs(tfgprev - tfgnew) > 2
%     cpfgnew = cpfg (tfgnew);
%     tfgprev = tfgnew;
%     tfgnew  = hf / cpfgnew;
% end
% tfg  = tfgnew                 % OK has been verified 
% cpfg = cpfg (tfg);

s0fg = EntropyFG(15);
sfg  = EntropyFG(tfg);
h0fg = EnthalpyFG(15) - EnthalpyFG(0);

exergyFG = Exergy(h0fg,hf,s0fg,sfg);

end

%%
function [ cpair ] = cpair (Ta)
% Function calculating the average specific heat for air between 0 and 15 °C

cpO20  = janaf('c','O2',273.15);
cpO2a  = janaf('c','O2',273.15 + Ta);
cpN20  = janaf('c','N2',273.15);
cpN2a  = janaf('c','N2',273.15 + Ta);

% cp for ambiant air
cpair = 0.5 * (0.21*(cpO20 + cpO2a) + 0.79*(cpN20 + cpN2a)); %kJ / kgair K

end

%%
function [ cpfg ] = cpfg (Tfg)
% Function calculating the average specific heat for flue gasses at
% temperature Tf in °C

Wfgasses  = (x*MWCO2 + 0.5*y*MWH2O + (lambda - 1)*MWO2 + (x + 0.25*(y-2*z))*3.76*MWN2); %kg            
        
CO2_frac = x*MWCO2/Wfgasses;             %kg_CO2 / kg_fg
H2O_frac = 0.5*y*MWH2O/Wfgasses;         %kg_H2O / kg_fg
O2_frac  = (lambda - 1)*MWO2/Wfgasses;
N2_frac  = (x + 0.25*(y-2*z))*3.76*MWN2/Wfgasses;

% Caracteristics of the components 
% 0 means at 0 °C, a means at ambiant temperature, f means at tfg

cpO20  = janaf('c','O2',273.15);
cpO2f  = janaf('c','O2',273.15 + Tfg);
cpN20  = janaf('c','N2',273.15);
cpN2f  = janaf('c','N2',273.15 + Tfg);
cpCO20 = janaf('c','CO2',273.15);
cpCO2f = janaf('c','CO2',273.15 + Tfg);
cpH2O0 = janaf('c','H2O',273.15);
cpH2Of = janaf('c','H2O',273.15 + Tfg);

% cp for flue gasses
cpfg = 0.5 * (CO2_frac*(cpCO20 + cpCO2f) + H2O_frac*(cpH2O0 + cpH2Of) ...
             + O2_frac*(cpO20  +  cpO2f) +  N2_frac*(cpN20  +  cpN2f)); %kJ / kgfg K
end

function [ tfg ] = TemperatureFG (hfg)
% Function calculating the temperature of the flue gasses

tfg0 = 1500;                  % First guess for T combustion
cpfg0 = cpfg (tfg0);          % First guess for cp flue gasses

tfgprev = tfg0;
tfgnew  = hfg / cpfg0;
while abs(tfgprev - tfgnew) > 2
    cpfgnew = cpfg (tfgnew);
    tfgprev = tfgnew;
    tfgnew  = hfg / cpfgnew;
end
tfg  = tfgnew;               % OK has been verified 
    
end
%%
function [ S ] = EntropyFG (Tfg)
% Function calculating the entropy of the flue gasses at temperature Tfg in °C

Wfgasses  = (x*MWCO2 + 0.5*y*MWH2O + (lambda - 1)*MWO2 + (x + 0.25*(y-2*z))*3.76*MWN2); %kg            
        
CO2_frac = x*MWCO2/Wfgasses; %kg_CO2 / kg_fg
H2O_frac = 0.5*y*MWH2O/Wfgasses; %kg_H2O / kg_fg
O2_frac  = (lambda - 1)*MWO2/Wfgasses;
N2_frac  = (x + 0.25*(y-2*z))*3.76*MWN2/Wfgasses;

SO2f  = janaf('s','O2',273.15 + Tfg);
SN2f  = janaf('s','N2',273.15 + Tfg);
SCO2f = janaf('s','CO2',273.15 + Tfg);
SH2Of = janaf('s','H2O',273.15 + Tfg);

% Entropy for flue gasses at T = Tfg
S = CO2_frac*SCO2f + H2O_frac*SH2Of + O2_frac*SO2f +  N2_frac*SN2f; %kJ / kgfg

end

%%
function [ hfg ] = EnthalpyFG (Tfg)
% Function calculating the enthalpy of the flue gasses at temperature Tfg in °C
%
% Errors of estimation are due to the janaf file. eg: calling
% janaf('h','H2',273) give a non-zero result, which doesn't agree with the
% conventions. 

Wfgasses  = (x*MWCO2 + 0.5*y*MWH2O + (lambda - 1)*MWO2 + (x + 0.25*(y-2*z))*3.76*MWN2); %kg            
        
CO2_frac = x*MWCO2/Wfgasses; %kg_CO2 / kg_fg
H2O_frac = 0.5*y*MWH2O/Wfgasses; %kg_H2O / kg_fg
O2_frac  = (lambda - 1)*MWO2/Wfgasses;
N2_frac  = (x + 0.25*(y-2*z))*3.76*MWN2/Wfgasses;

HO2f  = janaf('h','O2',273.15 + Tfg);
HN2f  = janaf('h','N2',273.15 + Tfg);
HCO2f = janaf('h','CO2',273.15 + Tfg);
HH2Of = janaf('h','H2O',273.15 + Tfg);

hfg = CO2_frac*HCO2f + H2O_frac*HH2Of + O2_frac*HO2f +  N2_frac*HN2f; %kJ / kgfg

end

%%
function [ ha ] = EnthalpyAir (T_air)
% Function calculating the enthalpy of air at temperature T_air in °C
%
% Errors of estimation are due to the janaf file. 
% eg: calling janaf('h','H2',273) give a non-zero result, which doesn't 
% agree with the conventions.           

HO2  = janaf('h','O2',273.15 + T_air);
HN2  = janaf('h','N2',273.15 + T_air);


ha = (0.21*MWO2*HO2 + 0.79*MWN2*HN2)/MWair; %kJ / kgair

end

%%
function ex = Exergy (h0,h,s0,s)
    t0 = 273.15 + 15; %[K]
    ex = h - h0 - t0*(s - s0);
end
end

