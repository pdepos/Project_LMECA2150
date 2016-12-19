function [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution] = SCAnalysis( state, P_el, FH, RH, x, y, z, t_exh, Ta, lambda )
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

p_bache0 = 2;

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
ma1     = ((MWO2+3.76*MWN2)*(x + (y - 2*z)/4))/(12*x + y + 16*z);  %[kg_air/kg_fuel]
Wfgasses  = (x*MWCO2 + 0.5*y*MWH2O + (lambda - 1)*0.5*y*MWO2 + (x + 0.25*(y-2*z))*3.76*MWN2); %kg            
        
CO2_frac = x*MWCO2/Wfgasses;             %kg_CO2 / kg_fg
H2O_frac = 0.5*y*MWH2O/Wfgasses;         %kg_H2O / kg_fg
O2_frac  = (lambda - 1)*0.5*y*MWO2/Wfgasses;
N2_frac  = (x + 0.25*(y-2*z))*3.76*MWN2/Wfgasses;

H0_CO2   = -393400/MWCO2;
H0_H2O   = -241800/MWH2O;
H0FG     = CO2_frac*H0_CO2 + H2O_frac*H0_H2O;

switch RH
    case 'off' 
        r = 0;
        if FH == 0
            [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] ...
                = RankineHirnAnalysis(state, P_el, t_exh);
        elseif FH > 0
            [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] ...
                = BasicFeedHeatingAnalysis(state, P_el, FH, t_exh, r);
        else
            warning('Negative number of feedheaters not allowed')
            EnergyDistribution = 0;
            FlowRates = 0;
        end
    case 'on'
        r = 1;
        if FH == 0
            [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] ...
                = BasicReHeatingAnalysis(state, P_el, t_exh);
        elseif FH > 0
            [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] ...
                = ReHeatingAndFHAnalysis(state, P_el, FH, r);
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
function [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] = RankineHirnAnalysis(state, P_el, t_exh)

ind1   = 1;
ind2   = 2;
ind3   = 5;
ind4   = 6;

% Vapor Flow Rate
wen_mov = state{ind3}.h - state{ind4}.h;
wen_op  = state{ind2}.h - state{ind1}.h;
wen_mcy = wen_mov - wen_op;
Q1    = (state{ind3}.h - state{ind2}.h);
Q2    = (state{ind4}.h - state{ind1}.h);

wex_mov = state{ind3}.e - state{ind4}.e;
wex_op  = state{ind2}.e - state{ind1}.e;
wex_mcy = wex_mov - wex_op;
DExG   = (state{ind3}.e - state{ind2}.e);
DExC   = (state{ind4}.e - state{ind1}.e);

m_dot_v = P_el / (wen_mov*eta_mec - wen_op/eta_mec);

% Finding the enthalpy of the FG at combustion temperature and at exhaust.
cpair_Ta = cpair(0,Ta); 
ha_Ta    = cpair_Ta*Ta; 
hc_15    = 15*cp15;
hfg_comb = (LHV + hc_15 + lambda*ma1*ha_Ta)/(lambda*ma1 + 1); 
Tcomb    = TemperatureFG (hfg_comb);
hfg_exh  = cpfg(0,t_exh)*t_exh;

epsilon_exh = ((lambda*ma1 + 1)*hfg_exh - lambda*ma1*ha_Ta)/LHV;
eta_gen     = 1 - 0.01 - epsilon_exh;
eta_cyclen  = wen_mcy / Q1;

m_dot_c = m_dot_v * Q1 / (eta_gen * LHV);
m_dot_f = (1 + lambda*ma1)*m_dot_c;
 
eta_toten  = P_el/(m_dot_c*LHV);

% Energetic Losses
Q_gen        = m_dot_v * Q1;     % Power recieved by the fluid at generator
Q_comb       = Q_gen / eta_gen;  % Primary Power needed from the fuel
ConEnLosses  = m_dot_v * Q2;
MecEnLosses  = m_dot_v*((1 - eta_mec) * wen_mov + (1/eta_mec - 1)*wen_op);
GenEnLosses  = (1 - eta_gen) * Q_comb;

% Exergetic Efficiencies
efg_comb = ExergyFlueGasses(Tcomb);  % is ok
efg_exh  = ExergyFlueGasses(t_exh);  % is ok
er       = ExergyR(ha_Ta,Ta);        % Exergy Reagents = 0 at T_reagents = 15°C


eta_combex  = m_dot_f * (efg_comb - er) / (m_dot_c * ec);
eta_chimnex = m_dot_f * (efg_comb - efg_exh) /(m_dot_f*(efg_comb - er));
eta_transex = m_dot_v * DExG /(m_dot_f*(efg_comb - efg_exh));
eta_gex     = eta_combex*eta_chimnex*eta_transex;
eta_rotex   = wen_mcy / wex_mcy;
eta_cyclex  = wen_mcy / DExG;
eta_totex   = P_el/(m_dot_c*ec);


% Exergetic Losses
MecExLosses = MecEnLosses;
ConExLosses = m_dot_v * DExC;
RotExLosses = (1 -   eta_rotex) * m_dot_v * wex_mcy;
TraExLosses = (1 - eta_transex) * m_dot_f * (efg_comb - efg_exh);
StaExLosses = (1 - eta_chimnex) * m_dot_f * (efg_comb - er);
ComExLosses = (1 -  eta_combex) * m_dot_c * ec;
          
% Outputs
FlowRates = {'m°v' 'm°c' 'm°f'; m_dot_v m_dot_c m_dot_f};
eta_en    = {'eta_mec', 'eta_gen', 'eta_cyclen', 'eta_toten';...
            eta_mec, eta_gen, eta_cyclen, eta_toten};
EnergyDistribution = ...
    {'Effective Power', 'Mechanical Losses', 'Condensor Losses', 'Generator Losses', 'Primary Energy';...
    P_el, MecEnLosses, ConEnLosses, GenEnLosses, Q_comb};
eta_ex = ...
    {'eta_mec' 'eta_totex' 'eta_gex' 'eta_combex' 'eta_chimnex' 'eta_transex' 'eta_rotex' 'eta_cyclex';...
    eta_mec eta_totex eta_gex eta_combex eta_chimnex eta_transex eta_rotex eta_cyclex};
ExergyDistribution = ...
    {'Effective Power' 'Mechanical Losses' 'Condensor Losses' 'Turbine & Pump Losses'...
    'Steam Generator Losses' 'Stack Losses' 'Combustion Losses';...
    P_el MecExLosses ConExLosses RotExLosses TraExLosses StaExLosses ComExLosses};
end

%%
function [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] = BasicFeedHeatingAnalysis(state, P_el, fh, t_exh, r)
  
alpha  = 10;
beta   = 4;  
gamma  = 4;
ind1   = 1;
ind2   = 2;
ind3   = 5;
ind4   = 6;
ind5   = 7;
ind6   = 8;
ind8   = 10;

%% Finding the Vapor Flow Rates
X = VFractions (state, fh, r);

if fh > gamma
    position_bache = 1; %References the roman number of the schematic
    index_bleed_bache = alpha + 1;
    while (state{index_bleed_bache}.p < p_bache0) && (position_bache <= fh)
        position_bache = position_bache + 1;
        index_bleed_bache = index_bleed_bache + beta;
    end
    mv = 1 + sum(X(1:position_bache - 1));
    ind6BP = alpha + beta*(position_bache + 1) - 2; % Index exit of the bache pump
    ind7BP = alpha + beta*position_bache - 1;       % Index entrance of bache pump
    wen_cp   = mv*abs(state{ind6}.h - state{ind5}.h);
    wen_bp   = (1 + sum(X))*abs(state{ind6BP}.h - state{ind7BP}.h);
    wen_fwp  = (1 + sum(X))*abs(state{ind2}.h - state{ind1}.h);
    wex_cp   = mv*abs(state{ind6}.e - state{ind5}.e);
    wex_bp   = (1 + sum(X))*abs(state{ind6BP}.e - state{ind7BP}.e);
    wex_fwp  = (1 + sum(X))*abs(state{ind2}.e - state{ind1}.e);
    Q2x      = mv*(state{ind8}.h - state{ind5}.h);
else
    mv       = (1 + sum(X));
    wen_cp   = mv*abs(state{ind6}.h - state{ind5}.h);
    wen_bp   = 0;
    wen_fwp  = mv*abs(state{ind2}.h - state{ind1}.h);
    wex_cp   = mv*abs(state{ind6}.e - state{ind5}.e);
    wex_bp   = 0;
    wex_fwp  = mv*abs(state{ind2}.e - state{ind1}.e);
    Q2x      = mv*abs(state{ind8}.h - state{ind5}.h);
end

wen_mov = (state{ind3}.h - state{ind4}.h);
for i = 1:fh
    ind4i = alpha + beta*(i-1) + 1; 
    wen_mov = wen_mov + X(i)*(state{ind3}.h - state{ind4i}.h);  
end
wen_op  = wen_cp + wen_bp + wen_fwp;
wen_mcy = wen_mov - wen_op;

Q1 = (1 + sum(X))*(state{ind3}.h - state{ind2}.h);
Q2 = (state{ind4}.h - state{ind5}.h) + Q2x;

wex_mcy = (state{ind3}.e - state{ind4}.e) - (wex_cp + wex_bp + wex_fwp);
for i = 1:fh
    ind4i = alpha + beta*(i-1) + 1; 
    wex_mcy = wex_mcy + X(i)*(state{ind3}.e - state{ind4i}.e);   
end
DExG = ((1 + sum(X))*(state{ind3}.e - state{ind2}.e));
DExC = ((state{ind4}.e - state{ind5}.e) + (mv - 1)*(state{ind8}.e - state{ind5}.e));

m_dot_vc = P_el /(wen_mov*eta_mec - wen_op/eta_mec);    % Vapor flow at the condensor

% Vapor Flow Rates
FlowRates      = FRCreation(fh + 4);
FlowRates{2,1} = m_dot_vc;
m_dot_vg       = m_dot_vc * (1 + sum(X));
FlowRates{2,2} = m_dot_vg;
for i = 1:fh
    index = 2 + i;
    FlowRates{2,index} = m_dot_vc*X(i);
end    

%% Energy Analysis and Other Flow Rates
% States of the Combustion 
if state{ind2}.t < t_exh
    tfg_exitvp = t_exh;
else
    tfg_exitvp = state{ind2}.t + 10;   % FG temperature after the water heat exchanger part of the boiler
end  
hc_15      = cp15 * 15;                % Fixing the fuel's temperature at 15 °C Its influence is very small in any case. 
hfg_exh    = cpfg(0,t_exh) * t_exh;
hfg_exitvp = cpfg(0,tfg_exitvp) * tfg_exitvp;
ha_ta      = cpair(0,Ta) * Ta;
ha_warm    = ha_ta + (1 + lambda*ma1)/(lambda*ma1) * (hfg_exitvp - hfg_exh);
hfg_comb   = (LHV*0.99 + hc_15 + lambda*ma1*ha_warm) / (lambda*ma1 + 1); % LHV*0.99 to account for 1% wall losses
Tcomb      = TemperatureFG(hfg_comb);

% epsilon_exh = ((lambda*ma1 + 1)*hfg_exh - lambda*ma1*ha_Ta)/LHV;
% eta_gen     = 1 - 0.01 - epsilon_exh; 
% m_dot_c  = m_dot_vc * Q1 / (eta_gen * LHV);
% m_dot_fg = (1 + lambda*ma1) * m_dot_c;
% FlowRates{2,fh + 4} = m_dot_fg;
% FlowRates{2,fh + 3} = m_dot_c;


m_dot_f = m_dot_vc * Q1 / (hfg_comb - hfg_exitvp);
m_dot_c  = m_dot_f / (1 + lambda*ma1);
FlowRates{2,fh + 4} = m_dot_f;
FlowRates{2,fh + 3} = m_dot_c;

eta_gen     = (hfg_comb - hfg_exh) / hfg_comb;
eta_cyclen  = wen_mcy / Q1;
eta_toten   = P_el / (m_dot_c*LHV);

Q_comb       = m_dot_c  * LHV;  % Primary Power delivered by the fuel
ConEnLosses  = m_dot_vc * Q2;
GenEnLosses  = (1 - eta_gen) * Q_comb;
MecEnLosses  = m_dot_vc * ((1 - eta_mec) * wen_mov + (1/eta_mec - 1) * wen_op);


%% Exergy Analysis
% Cycle Efficiency
eta_rotex  = wen_mcy / wex_mcy;
eta_cyclex = wen_mcy / DExG;

% Combustion Efficiency
efg_comb   = ExergyFlueGasses(Tcomb); 
efg_exh    = ExergyFlueGasses(t_exh);
efg_exitvp = ExergyFlueGasses(tfg_exitvp);  % Exergy of FG after heat transfer to vapor 
efg_ta      = ExergyFlueGasses(Ta);          % Exergy of the FG at ambiant atmosphere
ea_warm    = ExergyR(ha_warm,tfg_exitvp);   % Exergy Reagents after warmup. We neglige 
                                            % the fuel preheating, so it is basically the exergy of warm air
ea_ta      = ExergyR(ha_ta,Ta);             % Exergy of the air at ambiant conditions

eta_totex   = P_el/(m_dot_c*ec);
eta_combex  = m_dot_f*(efg_comb - ea_warm) / (m_dot_c * ec);
eta_chimnex = (efg_comb - efg_exh) / (efg_comb - efg_ta);
eta_transex = (m_dot_vc*DExG + m_dot_c*lambda*ma1*(ea_warm - ea_ta)) / (m_dot_f*(efg_comb - efg_exh));
eta_gex     = eta_combex*eta_chimnex*eta_transex;

% Exergetic Losses
MecExLosses  = MecEnLosses;
ConExLosses  = m_dot_vc * DExC;
ComExLosses  = (1 - eta_combex) * m_dot_c  * ec;
RotExLosses  = (1 -  eta_rotex) * m_dot_vc * wex_mcy ;
StaExLosses  = m_dot_f*efg_exh;
FGVTransExL  = m_dot_f*(efg_comb - efg_exitvp) - m_dot_vc * DExG;
FGATransExL  = m_dot_f*(efg_exitvp - efg_exh) - m_dot_c*lambda*ma1*(ea_warm - ea_ta);

% Heat Exchanger (HE) and valve losses
if fh <= gamma
    HeatExLosses = - m_dot_vg*(state{ind1}.e - state{ind6}.e);
    for i = 1:fh
        indFR = i + 2;
        ind4i = alpha + beta*(i - 1) + 1; 
        HeatExLosses = HeatExLosses + FlowRates{2,indFR}*(state{ind4i}.e - state{ind8}.e);
    end
elseif fh > gamma
    HeatExLosses = - m_dot_vg * (state{ind1}.e - state{ind6BP}.e)...
                   - m_dot_vc * mv * (state{ind7BP}.e - state{ind6}.e);
    for i = 1:fh
    indFR     = i + 2;
    ind4i     = alpha + beta*(i - 1) + 1;     
        if i < position_bache
            HeatExLosses = HeatExLosses + FlowRates{2,indFR}*(state{ind4i}.e - state{ind8}.e);
        elseif i >= position_bache
            HeatExLosses = HeatExLosses + FlowRates{2,indFR}*(state{ind4i}.e - state{ind7BP}.e);
        end
    end
end

%% Outputs
eta_en    = {'eta_mec', 'eta_gen', 'eta_cyclen', 'eta_toten';...
              eta_mec, eta_gen, eta_cyclen, eta_toten};
EnergyDistribution = ...
    {'Effective Power', 'Mechanical Losses', 'Condensor Losses', 'Generator Losses', 'Primary Energy';...
    P_el, MecEnLosses, ConEnLosses, GenEnLosses, Q_comb};
eta_ex = ...
    {'eta_mec' 'eta_totex' 'eta_gex' 'eta_combex' 'eta_chimnex' 'eta_transex' 'eta_rotex' 'eta_cyclex';...
    eta_mec eta_totex eta_gex eta_combex eta_chimnex eta_transex eta_rotex eta_cyclex};
ExergyDistribution = ...
    {'Effective Power' 'Mechanical Losses' 'Condensor Losses' 'Turbine & Pump Losses'...
     'Steam Generator Losses' 'Air Preheating Losses' 'FW Heater Losses' 'Stack Losses' 'Combustion Losses'; ...
     P_el MecExLosses ConExLosses RotExLosses FGVTransExL FGATransExL HeatExLosses StaExLosses ComExLosses};

end

%%
function [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] = BasicReHeatingAnalysis(state, P_el, t_exh)
% Function Analysing the BasicReHeating Cycle.

ind1    = 1;
ind2    = 2;
ind3HP  = 5;
ind4HP  = 6;
ind3LP  = 7;
ind4LP  = 8;


wen_mov = (state{ind3HP}.h - state{ind4HP}.h) + (state{ind3LP}.h - state{ind4LP}.h);
wen_op  = (state{ind2}.h - state{ind1}.h);
wen_mcy = wen_mov - wen_op;
Q1      = (state{ind3HP}.h - state{ind2}.h) + (state{ind3LP}.h - state{ind4HP}.h);
Q2      = (state{ind4LP}.h - state{ind1}.h);

wex_mcy = (state{ind3HP}.e - state{ind4HP}.e) + (state{ind3LP}.e - state{ind4LP}.e) - ...
          (state{ind2}.e - state{ind1}.e);
DExG    = (state{ind3HP}.e - state{ind2}.e) + (state{ind3LP}.e - state{ind4HP}.e);
DExC    = (state{ind4LP}.e - state{ind1}.e);

%% Vapor Flow Rates
m_dot_v = P_el / (wen_mov*eta_mec - wen_op/eta_mec);

%% Energy Analysis
% Finding the enthalpy of the FG at combustion temperature and at exhaust.
hc_15    = 15*cp15;
ha_Ta    = cpair(0,Ta)*Ta; 
hfg_exh  = cpfg(0,t_exh)*t_exh;
hfg_comb = (LHV*0.99 + hc_15 + lambda*ma1*ha_Ta)/(lambda*ma1 + 1);
Tcomb    = TemperatureFG (hfg_comb); 

% Efficiencies and Combustion Flow Rates
epsilon_exh = ((lambda*ma1 + 1)*hfg_exh - lambda*ma1*ha_Ta)/LHV;
eta_gen     = 1 - 0.01 - epsilon_exh; 
eta_cyclen  = wen_mcy / Q1;

m_dot_c  = m_dot_v * Q1 / (eta_gen * LHV);
m_dot_f = (1 + lambda*ma1)*m_dot_c; 

eta_toten  = P_el / (m_dot_c*LHV);

% Losses
Q_gen        = m_dot_v * Q1;     % Power recieved by the fluid at generator
Q_comb       = Q_gen / eta_gen;  % Primary Power needed from the fuel
ConEnLosses  = m_dot_v * Q2;
MecEnLosses  = m_dot_v*((1 - eta_mec) * wen_mov + (1/eta_mec - 1)*wen_op);%(1 - eta_mec) * P_el;
GenEnLosses  = (1 - eta_gen) * Q_comb;

% Exergetic Efficiencies
efg_comb = ExergyFlueGasses(Tcomb);   % is ok
efg_exh  = ExergyFlueGasses(t_exh);   % is ok
er       = ExergyR(ha_Ta,Ta);         % Exergy Reagents = 0 at T_reagents = 15°C

eta_totex   = P_el/(m_dot_c*ec);
eta_combex  = m_dot_f * (efg_comb - er) / (m_dot_c * ec);
eta_chimnex = m_dot_f * (efg_comb - efg_exh) /(m_dot_f*(efg_comb - er));
eta_transex = m_dot_v * DExG /(m_dot_f*(efg_comb - efg_exh));
eta_gex     = eta_combex*eta_chimnex*eta_transex;
eta_rotex   = wen_mcy / wex_mcy;
eta_cyclex  = wen_mcy / DExG;
          
% Exergetic Losses
MecExLosses = MecEnLosses;
ConExLosses = m_dot_v * DExC;
RotExLosses = (1 - eta_rotex)   * m_dot_v *  wex_mcy;
TraExLosses = (1 - eta_transex) * m_dot_f * (efg_comb - efg_exh);
StaExLosses = (1 - eta_chimnex) * m_dot_f * (efg_comb - er);
ComExLosses = (1 -  eta_combex) * m_dot_c * ec;
          
% Outputs
FlowRates = {'m°v' 'm°c' 'm°f'; m_dot_v m_dot_c m_dot_f};
eta_en    = {'eta_mec', 'eta_gen', 'eta_cyclen', 'eta_toten';...
            eta_mec, eta_gen, eta_cyclen, eta_toten};
EnergyDistribution = ...
    {'Effective Power', 'Mechanical Losses', 'Condensor Losses', 'Generator Losses', 'Primary Energy';...
    P_el, MecEnLosses, ConEnLosses, GenEnLosses, Q_comb};
eta_ex = ...
    {'eta_mec' 'eta_totex' 'eta_gex' 'eta_combex' 'eta_chimnex' 'eta_transex' 'eta_rotex' 'eta_cyclex';...
    eta_mec eta_totex eta_gex eta_combex eta_chimnex eta_transex eta_rotex eta_cyclex};
ExergyDistribution = ...
    {'Effective Power' 'Mechanical Losses' 'Condensor Losses' 'Turbine & Pump Losses'...
    'Steam Generator Losses' 'Stack Losses' 'Combustion Losses';...
    P_el MecExLosses ConExLosses RotExLosses TraExLosses StaExLosses ComExLosses};
end

%%
function [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] = ReHeatingAndFHAnalysis(state, P_el, fh, r)
alpha = 12;                %Number of fixed points of the schematic
beta  = 4;                 %Number of points per state: eg: beta = 3 means 4.1 6.1 7.1 or 4.2 6.2 7.2
gamma = 4;
ind1    = 1;
ind2    = 2;
ind3HP  = 5;
ind4HP  = 6;
ind3LP  = 7;
ind4LP  = 8;
ind5    = 9;
ind6    = 10;
ind8    = 12;

%% Finding the Vapor Flow Rates
X = VFractions (state, fh, r);

if fh > gamma
    position_bache = 1;       %References the roman number of the schematic
    index_bleed_bache = alpha + 1;
    while (state{index_bleed_bache}.p < p_bache0) && (position_bache <= fh)
        position_bache = position_bache + 1;
        index_bleed_bache = index_bleed_bache + beta;
    end
    mv = 1 + sum(X(1:position_bache - 1));
    ind6BP = alpha + beta*(position_bache + 1) - 2; % Index exit of the bache pump
    ind7BP = alpha + beta*position_bache - 1;       % Index entrance of bache pump
    wen_cp   = mv*abs(state{ind6}.h - state{ind5}.h);
    wen_bp   = (1 + sum(X))*abs(state{ind6BP}.h - state{ind7BP}.h);
    wen_fwp  = (1 + sum(X))*abs(state{ind2}.h - state{ind1}.h);
    wex_cp   = mv*abs(state{ind6}.e - state{ind5}.e);
    wex_bp   = (1 + sum(X))*abs(state{ind6BP}.e - state{ind7BP}.e);
    wex_fwp  = (1 + sum(X))*abs(state{ind2}.e - state{ind1}.e);
    Q2x      = mv*(state{ind8}.h - state{ind5}.h);
else
    mv       = (1 + sum(X));
    wen_cp   = mv*abs(state{ind6}.h - state{ind5}.h);
    wen_bp   = 0;
    wen_fwp  = mv*abs(state{ind2}.h - state{ind1}.h);
    wex_cp   = mv*abs(state{ind6}.e - state{ind5}.e);
    wex_bp   = 0;
    wex_fwp  = mv*abs(state{ind2}.e - state{ind1}.e);
    Q2x      = mv*abs(state{ind8}.h - state{ind5}.h);
end

wen_mov = (1 + sum(X))*(state{ind3HP}.h - state{ind4HP}.h) + (state{ind3LP}.h - state{ind4LP}.h);
for i = 1:fh
    ind4i = alpha + beta*(i-1) + 1; 
    wen_mov = wen_mov + X(i)*(state{ind3LP}.h - state{ind4i}.h);  
end
wen_op = wen_cp + wen_bp + wen_fwp;
wen_mcy = wen_mov - wen_op;

Q1 = (1 + sum(X))*(state{ind3HP}.h - state{ind2}.h) + (1 + sum(X(1:fh)))*(state{ind3LP}.h - state{ind4HP}.h);
Q2 = (state{ind4LP}.h - state{ind5}.h) + Q2x;

wex_mcy = (1 + sum(X))*(state{ind3HP}.e - state{ind4HP}.e) + (state{ind3LP}.e - state{ind4LP}.e)...
        - wex_cp - wex_bp - wex_fwp;
for i = 1:fh
    ind4i = alpha + beta*(i-1) + 1; 
    wex_mcy = wex_mcy + X(i)*(state{ind3LP}.e - state{ind4i}.e);   
end
DExG = (1 + sum(X))*(state{ind3HP}.e - state{ind2}.e) + (1 + sum(X(1:fh)))*(state{ind3LP}.e - state{ind4HP}.e);
DExC = (state{ind4LP}.e - state{ind5}.e + (mv - 1)*(state{ind8}.e - state{ind5}.e));

m_dot_vc = P_el / (wen_mov*eta_mec - wen_op/eta_mec);   % Vapor flow in condensor

% Vapor Flow Rates
FlowRates      = FRCreation(fh + 5);
FlowRates{2,1} = m_dot_vc;
m_dot_vg       = m_dot_vc * (1 + sum(X));
FlowRates{2,2} = m_dot_vg;
for i = 1:(fh + 1)
    index = 2 + i;
    FlowRates{2,index} = FlowRates{2,1}*X(i);
end    

%% Energy Analysis and Other Flow Rates
% States of the Combustion 
if state{ind2}.t < t_exh
    tfg_exitvp = t_exh;
else
    tfg_exitvp = state{ind2}.t + 10;   % FG temperature after the water heat exchanger part of the boiler
end  
hc_15      = 15*cp15;                  % Fixing the fuel's temperature at 15 °C Its influence is very small in any case. 
hfg_exh    = cpfg(0,t_exh)*t_exh;
hfg_exitvp = cpfg(0,tfg_exitvp)*tfg_exitvp;
ha_ta      = cpair(0,Ta)*Ta;
ha_warm    = ha_ta + (1 + lambda*ma1)/(lambda*ma1)*(hfg_exitvp - hfg_exh);
hfg_comb   = (LHV*0.99 + hc_15 + lambda*ma1*ha_warm)/(lambda*ma1 + 1);      % LHV*0.99 to account for 1% wall losses
Tcomb      = TemperatureFG (hfg_comb);

m_dot_f  = m_dot_vc * Q1 / (hfg_comb - hfg_exitvp);
m_dot_c  = m_dot_f / (1 + lambda*ma1);
FlowRates{2,fh + 5} = m_dot_f;
FlowRates{2,fh + 4} = m_dot_c;

eta_cyclen = wen_mcy / Q1;
eta_toten  = P_el / (m_dot_c*LHV);
eta_gen    = (hfg_comb - hfg_exh) / hfg_comb;

Q_comb      = m_dot_c  * LHV;  % Primary Power delivered by the fuel
ConEnLosses = m_dot_vc * Q2;
GenEnLosses = (1 - eta_gen) * Q_comb;
MecEnLosses = m_dot_vc * ((1 - eta_mec) * wen_mov + (1/eta_mec - 1) * wen_op);

%% Exergy Analysis
% Cycle Efficiency
eta_rotex  = wen_mcy / wex_mcy;
eta_cyclex = wen_mcy / DExG;

% Combustion Efficiency
efg_comb   = ExergyFlueGasses(Tcomb); 
efg_exh    = ExergyFlueGasses(t_exh);
efg_exitvp = ExergyFlueGasses(tfg_exitvp);  % Exergy of FG after heat transfer to vapor 
efg_ta     = ExergyFlueGasses(Ta);          % Exergy of the FG at ambiant atmosphere
ea_warm    = ExergyR(ha_warm,tfg_exitvp);   % Exergy Reagents after warmup. We neglige 
                                            % the fuel preheating, so it is basically the exergy of warm air
ea_ta      = ExergyR(ha_ta,Ta);             % Exergy of the air at ambiant conditions

eta_totex   = P_el / (m_dot_c*ec);
eta_combex  = m_dot_f * (efg_comb - ea_warm) / (m_dot_c * ec);
eta_chimnex = (efg_comb - efg_exh) / (efg_comb - efg_ta);
eta_transex = (m_dot_vc*DExG + m_dot_c*lambda*ma1*(ea_warm - ea_ta)) / (m_dot_f*(efg_comb - efg_exh));
eta_gex     = eta_combex*eta_chimnex*eta_transex;


% Exergetic Losses
MecExLosses  = MecEnLosses;
ConExLosses  = m_dot_vc * DExC;
ComExLosses  = (1 - eta_combex) * m_dot_c  * ec;
RotExLosses  = (1 -  eta_rotex) * m_dot_vc * wex_mcy ;
StaExLosses  = m_dot_f*efg_exh;
FGVTransExL  = m_dot_f*(efg_comb - efg_exitvp) - m_dot_vc * DExG;
FGATransExL  = m_dot_f*(efg_exitvp - efg_exh) - m_dot_c*lambda*ma1*(ea_warm - ea_ta);

% Heat Exchanger (HE) losses and valve losses
if fh <= gamma
    HeatExLosses = - m_dot_vg*(state{ind1}.e - state{ind6}.e);
    for i = 1:fh + 1
        indFR = i + 2;
        ind4i = alpha + beta*(i - 1) + 1; 
        HeatExLosses = HeatExLosses + FlowRates{2,indFR}*(state{ind4i}.e - state{ind8}.e);
    end
elseif fh > gamma
    HeatExLosses = - m_dot_vg * (state{ind1}.e - state{ind6BP}.e)...
                   - m_dot_vc * mv * (state{ind7BP}.e - state{ind6}.e);
    for i = 1:fh + 1
    indFR     = i + 2;
    ind4i     = alpha + beta*(i - 1) + 1;     
        if i < position_bache
            HeatExLosses = HeatExLosses + FlowRates{2,indFR}*(state{ind4i}.e - state{ind8}.e);
        elseif i >= position_bache
            HeatExLosses = HeatExLosses + FlowRates{2,indFR}*(state{ind4i}.e - state{ind7BP}.e);
        end
    end
end

%% Outputs
eta_en    = {'eta_mec', 'eta_gen', 'eta_cyclen', 'eta_toten';...
              eta_mec, eta_gen, eta_cyclen, eta_toten};
EnergyDistribution = ...
    {'Effective Power', 'Mechanical Losses', 'Condensor Losses', 'Generator Losses', 'Primary Energy';...
    P_el, MecEnLosses, ConEnLosses, GenEnLosses, Q_comb};
eta_ex = ...
    {'eta_mec' 'eta_totex' 'eta_gex' 'eta_combex' 'eta_chimnex' 'eta_transex' 'eta_rotex' 'eta_cyclex';...
    eta_mec eta_totex eta_gex eta_combex eta_chimnex eta_transex eta_rotex eta_cyclex};
ExergyDistribution = ...
    {'Effective Power' 'Mechanical Losses' 'Condensor Losses' 'Turbine & Pump Losses'...
     'Steam Generator Losses' 'Air Preheating Losses' 'FW Heater Losses' 'Stack Losses' 'Combustion Losses'; ...
     P_el MecExLosses ConExLosses RotExLosses FGVTransExL FGATransExL HeatExLosses StaExLosses ComExLosses};

end

%%
function [ exergyFG ] = ExergyFlueGasses(tfg)
% Function Calculating the Exergy of the Flue Gasses temperature tfg in [°C]. 

% Old Version

delta_sfg  = cplogfg(15,tfg)*log((tfg + 273.15)/(15 + 273.15));
delta_hfg  = cpfg(15,tfg)*(tfg - 15);

t0 = 273.15 + 15; %[K]
exergyFG = delta_hfg - t0*delta_sfg;
end

%% 
function [ er ] = ExergyR (ha_warm,t)
% Fonction calculating the exergy of the reagents with 't' the temperature of 
% the reagents in °C or a first guess temperature to solve for ha_warm.

Twarm0 = t;
cpa0   = cpair(0,Twarm0);

i = 0;
Twarm  = ha_warm/cpa0;
error = Twarm - Twarm0;
while (abs(error) > 1) && (i < 50)
    cpa = cpair(0,Twarm);
    Twarm0 = Twarm;
    Twarm  = ha_warm/cpa;
    error  = Twarm - Twarm0;
    i = i + 1;
    if i >= 49
        warning ('max number of iterations in ExergyR reached')
    end
end

t0 = 273.15 + 15; 
ea = cpair(15,Twarm)*(Twarm - 15) - t0*cplogair(Twarm,15)*log((Twarm + 273.15)/(15 + 273.15));
er = ea * lambda*ma1/(1 + lambda*ma1);
end

%%
function [ cpair ] = cpair (t1,t2)
% Function calculating the average specific heat for air between t1 and t2 in °C
% Addapted the amount of N2 in the air composition to have a more precise value
% for cp, cross checked with table values.

T = linspace(273.15 + t1,273.15 + t2);
cpair = mean(0.21*janaf('c','O2',T) + 0.782*janaf('c','N2',T));

end

%%
function cplog = cplogair(t1, t2)
% Function calculating the logarythmic cp of air between t1 and t2 in °C
t1 = t1 + 273.15;
t2 = t2 + 273.15;

if t1 == t2
    cplog = 0.21*janaf('c','O2',t1) + 0.782*janaf('c','N2',t1);
else
    fO2 = (@(T) janaf('c','O2',T)./T);  
    fN2 = (@(T) janaf('c','N2',T)./T);

    cpO2 = integral(fO2 ,t1 ,t2)/log(t2/t1); 
    cpN2 = integral(fN2 ,t1 ,t2)/log(t2/t1);  

    cplog = 0.21*cpO2 + 0.782*cpN2;
end
end

%%
function [ cpfg ] = cpfg (t1,t2)
% Function calculating the average specific heat for flue gasses between 
% temperatures T1 and T2 in [°C]

T = linspace(273.15 + t1,273.15 + t2);
cpfg = mean(CO2_frac*janaf('c','CO2',T) + H2O_frac*janaf('c','H2O',T)...
         + O2_frac*janaf('c','O2',T)  +  N2_frac*janaf('c','N2',T));
end

%%
function [ cplog ] = cplogfg (t1, t2)
% Function calculating the logarythmic cp of the flue gasses between t1 and t2 in °C

t1 = t1 + 273.15;
t2 = t2 + 273.15;
if t1 == t2
    cplog =  O2_frac*janaf('c','O2' ,t1) +  N2_frac*janaf('c','N2' ,t1) +...
            CO2_frac*janaf('c','CO2',t1) + H2O_frac*janaf('c','H2O',t1);
else
    fO2  = (@(T) janaf('c','O2' ,T)./T);  
    fN2  = (@(T) janaf('c','N2' ,T)./T);
    fCO2 = (@(T) janaf('c','CO2',T)./T);
    fH2O = (@(T) janaf('c','H2O',T)./T);


    cpO2  = integral(fO2  ,t1 ,t2)/log(t2/t1); 
    cpN2  = integral(fN2  ,t1 ,t2)/log(t2/t1);  
    cpCO2 = integral(fCO2 ,t1 ,t2)/log(t2/t1); 
    cpH2O = integral(fH2O ,t1 ,t2)/log(t2/t1); 

    cplog = O2_frac*cpO2 + N2_frac*cpN2 + CO2_frac*cpCO2 + H2O_frac*cpH2O;
end

end

%%
function [ tfg ] = TemperatureFG (hfg)
% Function calculating the adiabatic flame temperature of the flue gasses. 
% Take the combustion enthalpy as parameter in [kJ/kg]

tfg0 = 1500;                % First guess for T combustion
cpfg0 = cpfg (0,tfg0);      % First guess for cp flue gasses

i = 0;
tfgnew = hfg / cpfg0;
error  = tfg0 - tfgnew;
while (abs(error) > 1) && (i < 50)
    cpfg = cpfg (0,tfgnew);
    tfg0 = tfgnew;
    tfgnew = hfg / cpfg;
    error  = tfg0 - tfgnew;
    i = i + 1;
    if i >= 49
        warning ('max number of iterations in TemperatureFG reached')
    end
end
tfg  = tfgnew;        
    
end

%%
function [ flowrates ] = FRCreation(fr)
% Function initiating the flowrate structure, with the number of flow rates
% as input argument. 

flowrates = cell(2,fr);
flowrates{1,1}    = 'm°vc';
flowrates{1,2}    = 'm°vg';
for i = 3:fr-2
    formatSpec1 = '%s%d';
    str1 = sprintf(formatSpec1,'m°vb',i-2);
    flowrates{1,i} = str1;    
end
flowrates{1,fr-1} = 'm°c';
flowrates{1,fr}   = 'm°f';

for i = 1:fr
    flowrates{2,i} = 0;
end
end

%%
function [X] = VFractions (state, fh, rh_param) 
% Function calculating the vapor fractions of cycles with feedheating
% enabled. Input parameter are:
% - The states of the cycle
% - The number of feedheaters
% - A parameter that notifies if there's Reheating or not. 

% Verifying if RH is on and adapting the number of bleeds
if rh_param == 0
    alpha = 10;
    beta  = 4;  
    gamma = 4;
    ind1   = 1;
    ind6   = 8;
    ind7   = 9;
    ind41  = 11;
elseif rh_param == 1
    fh = fh + 1;
    alpha = 12;               
    beta  = 4; 
    gamma = 5;
    ind1    = 1;
    ind6    = 10;
    ind7    = 11;
    ind41   = 13;
else
    warning ('Problems in VFractions with rh_param')
end
%% Finding the Vapor Flow Rates

% Left Hand Side of the flow rates system
A = zeros(fh,1);
% Right Hand Side of the flow rates system
B = zeros(fh);

position_bache = 1;                   %References the roman number of the schematic
if fh == 1
    A(1) = state{ind1}.h - state{ind6}.h;
    B(1) = state{ind41}.h - state{ind7}.h - A(1);
elseif fh > gamma                     %Finding the bache if the number of bleeds is higher than gamma
    index_bleed_bache = alpha + 1;
    while (state{index_bleed_bache}.p < p_bache0) && (position_bache <= fh)
        position_bache = position_bache + 1;
        index_bleed_bache = index_bleed_bache + beta;
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
            elseif (i >  position_bache) && (j > i)
                B(i,j) = (state{ind8iplus}.h - state{ind7i}.h) - A(i);
            
            % To detect errors
            else
                warning ('Error during VFraction Matrix B coefficient calculation')
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
        else
            A(i) = state{ind6iplus}.h - state{ind6i}.h;
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
                warning ('Error during VFraction Matrix B coefficient calculation')
                B(i,j) = NaN;
            end
        end
    end
end

X = B\A;
end

end

