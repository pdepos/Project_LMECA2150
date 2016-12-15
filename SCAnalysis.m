function [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution] = SCAnalysis( state, P_el_LP, P_el_HP, FH, RH, x, y, z, t_exh, Ta, lambda )
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
function [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] = RankineHirnEnergy(state, P_el, t_exh)

ind1   = 1;
ind2   = 2;
ind3   = 5;
ind4   = 6;

% Vapor Flow Rate
m_dot_v = P_el / (((state{ind3}.h - state{ind4}.h) - (state{ind2}.h - state{ind1}.h))*eta_mec);

% Finding the enthalpy of the FG at combustion temperature and at exhaust.
cpair_Ta = cpair(0,Ta); 
ha_Ta    = cpair_Ta*Ta; 
hc_15    = 15*cp15;
hfg_comb = (LHV*0.99 + hc_15 + lambda*ma1*ha_Ta)/(lambda*ma1 + 1); % Multiply LHV by 0.99 to account for 1% wall losses
Tcomb    = TemperatureFG (hfg_comb);
hfg_exh  = cpfg(0,t_exh)*t_exh;

% With this, I find the flow rate of FG
m_dot_f = m_dot_v*(state{ind3}.h - state{ind2}.h)/(hfg_comb - hfg_exh);

% And flow rate of combustible
m_dot_c = m_dot_f/(1 + lambda*ma1);

% Energetic Efficiencies
eta_cyclen = ((state{ind3}.h - state{ind4}.h) - abs(state{ind2}.h - state{ind1}.h)) / ...
             (state{ind3}.h - state{ind2}.h);
eta_gen    =  m_dot_v*(state{ind3}.h - state{ind2}.h) / (m_dot_c*LHV);  
eta_toten  = eta_mec*eta_gen*eta_cyclen;

% Energetic Losses
Q_gen        = m_dot_v * (state{ind3}.h - state{ind2}.h); % Power recieved by the fluid at generator
Q_comb       = Q_gen / eta_gen;                           % Primary Power needed from the fuel
ConEnLosses  = m_dot_v * (state{ind4}.h - state{ind1}.h);
MecEnLosses  = (1 - eta_mec) * P_el;
GenEnLosses  = (1 - eta_gen) * Q_comb;

% Exergetic Efficiencies
efg_comb = ExergyFlueGasses(Tcomb);  % is ok
efg_exh  = ExergyFlueGasses(t_exh);   % is ok
er       = ExergyR(ha_Ta,Ta);                  % Exergy Reagents = 0 at T_reagents = 15°C

eta_totex   = P_el/(m_dot_c*ec);
eta_combex  = m_dot_f*(efg_comb - er) / (m_dot_c * ec);
eta_chimnex = m_dot_f*(efg_comb - efg_exh)/(m_dot_f*(efg_comb - er));
eta_transex = m_dot_v*(state{ind3}.e - state{ind2}.e)/(m_dot_f*(efg_comb - efg_exh));
eta_gex     = eta_combex*eta_chimnex*eta_transex;
%eta_gextest = m_dot_v*(state{ind3}.e - state{ind2}.e) / (m_dot_c * ec);
eta_rotex   = ((state{ind3}.h - state{ind4}.h) - (state{ind2}.h - state{ind1}.h)) /...
              ((state{ind3}.e - state{ind4}.e) - (state{ind2}.e - state{ind1}.e));
eta_cyclex  = ((state{ind3}.h - state{ind4}.h) - (state{ind2}.h - state{ind1}.h)) /...
              ( state{ind3}.e - state{ind2}.e);
          
% Exergetic Losses
MecExLosses = MecEnLosses;
ConExLosses = m_dot_v * (state{ind4}.e - state{ind1}.e);
RotExLosses = (1 -   eta_rotex) * m_dot_v * ((state{ind3}.e - state{ind4}.e) - (state{ind2}.e - state{ind1}.e));
TraExLosses = (1 - eta_transex) * m_dot_f * (efg_comb - efg_exh);
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
end

%%
function [ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] = BasicFeedHeatingEnergy(state, P_el, fh, t_exh, r)

    
alpha = 10;
beta  = 4;  
gamma = 4;

ind1   = 1;
ind2   = 2;
ind3   = 5;
ind4   = 6;
ind5   = 7;
ind6   = 8;
ind7   = 9;
ind8   = 10;
ind41  = 11;

n = alpha + beta * fh;

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
    while state{index_bleed_bache}.p < 4.6
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

X = B\A;
w_turbine = state{ind3}.h - state{ind4}.h;
for i = 1:fh
    ind4i     = alpha + beta*(i-1) + 1; 
    w_turbine = w_turbine + X(i)*(state{ind3}.h - state{ind4i}.h);  
end

P_meca = P_el/eta_mec;
m_dot_vc = P_meca / w_turbine;     % Vapor flow in condensor

% Vapor Flow Rates
FlowRates      = FRCreation(fh + 4);
FlowRates{2,1} = m_dot_vc;
FlowRates{2,2} = FlowRates{2,1};
for i = 1:fh
    index = 2 + i;
    FlowRates{2,index} = FlowRates{2,1}*X(i);
    FlowRates{2,2}     = FlowRates{2,2} + FlowRates{2,index};
end
m_dot_vg = FlowRates{2,2};        % Vapor flow in generator

%% Energy Analysis and Other Flow Rates

% States of the Combustion 
if state{ind1}.t < t_exh
    tfg_exitvp = t_exh;
else
    tfg_exitvp = state{ind1}.t + 10;   % FG temperature after the water heat exchanger part of the boiler
end  
hc_15      = 15*cp15;                  % Fixing the fuel's temperature at 15 °C Its influence is very small in any case. 
hfg_exh    = cpfg(0,t_exh)*t_exh;
hfg_exitvp = cpfg(0,tfg_exitvp)*tfg_exitvp;
ha_Ta      = cpair(0,Ta)*Ta;
ha_warm    = ha_Ta + (1 + lambda*ma1)/(lambda*ma1)*(hfg_exitvp - hfg_exh);
hfg_comb   = (LHV*0.99 + hc_15 + lambda*ma1*ha_warm)/(lambda*ma1 + 1); % Multiply LHV by 0.99 to account for 1% wall losses
Tcomb      = TemperatureFG (hfg_comb);

% FR flue gasses and fuel
m_dot_fg = FlowRates{2,2}*(state{ind3}.h - state{ind2}.h)/(hfg_comb - hfg_exitvp);
m_dot_c  = m_dot_fg /(1 + lambda*ma1);
FlowRates{2,fh + 4} = m_dot_fg;
FlowRates{2,fh + 3} = m_dot_c;

if fh > gamma
    mv = m_dot_vc;   % Amount of water flowing through condensor pump
    for i = 1:position_bache - 1
        ind = i + 2;
        mv  = mv + FlowRates{2,ind};
    end
    ind6BP = alpha + beta*(position_bache + 1) - 2; % Index exit of the bache pump
    ind7BP   = alpha + beta*position_bache - 1;       % Index entrance of bache pump
    P_cp      = mv*abs(state{ind6}.h - state{ind5}.h);
    P_bp      = m_dot_vg*abs(state{ind6BP}.h - state{ind7BP}.h);
    P_fwp     = m_dot_vg*abs(state{ind2}.h - state{ind1}.h);
else
    mv = m_dot_vg;
    P_cp      = mv*abs(state{ind6}.h - state{ind5}.h);
    P_bp      = 0;
    P_fwp     = mv*abs(state{ind2}.h - state{ind1}.h);
end

eta_cyclen = (P_el/eta_mec - P_cp - P_bp - P_fwp) / (m_dot_vg*(state{ind3}.h - state{ind2}.h));
eta_gen    = m_dot_vg*(state{ind3}.h - state{ind2}.h) / (m_dot_c*LHV);  
eta_toten  = eta_mec*eta_gen*eta_cyclen;

% Energetic Losses
Q_gen        = m_dot_vg * (state{ind3}.h - state{ind2}.h); % Power recieved by the fluid at generator
Q_comb       = Q_gen / eta_gen;                            % Primary Power needed from the fuel
ConEnLosses  = m_dot_vc * (state{ind4}.h - state{ind1}.h);
MecEnLosses  = (1 - eta_mec) * P_el;
GenEnLosses  = (1 - eta_gen) * Q_comb;


%% Exergy Analysis
% Cycle Efficiency
denrotex = m_dot_vc*(state{ind3}.e - state{ind4}.e);
for i = 1:fh
    ind   = i + 2;
    ind4i = alpha + beta*(i-1) + 1; 
    denrotex = denrotex + FlowRates{2,ind}*(state{ind3}.e - state{ind4i}.e);   
end
denrotex = denrotex - m_dot_vg*((state{ind6}.e - state{ind5}.e) + (state{ind2}.e - state{ind1}.e));
eta_rotex  = (P_el/eta_mec - P_cp - P_bp - P_fwp)/denrotex;
eta_cyclex = (P_el/eta_mec - P_cp - P_bp - P_fwp)/(m_dot_vg*(state{ind3}.e - state{ind2}.e));

% Combustion Efficiency
efg_comb   = ExergyFlueGasses(Tcomb); 
efg_exh    = ExergyFlueGasses(t_exh);
efg_exitvp = ExergyFlueGasses(tfg_exitvp);     % Exergy of FG after heat transfer to vapor 
er         = ExergyR(ha_warm,tfg_exitvp);      % Exergy Reagents after warmup. We do not pre-heat the fuel.
ea         = ExergyR(ha_Ta,Ta);                % Exergy of ambiant atmosphere

eta_combex  = m_dot_fg*(efg_comb - er) / (m_dot_c * ec);
eta_chimnex = m_dot_fg*(efg_comb - efg_exh)/(m_dot_fg*(efg_comb - ea));
eta_transex = m_dot_vg*(state{ind3}.e - state{ind2}.e)/(m_dot_fg*(efg_comb - efg_exitvp));
eta_gex     = eta_combex*eta_chimnex*eta_transex;

eta_totex   = eta_mec*eta_gex*eta_cyclex;

% Exergetic Losses
MecExLosses  = MecEnLosses;
ConExLosses  = m_dot_vc * (state{ind4}.e - state{ind5}.e) + (mv - m_dot_vc)*(state{ind8}.e - state{ind5}.e);
RotExLosses  = (1 -   eta_rotex) * denrotex;
TraExLosses  = (1 - eta_transex) * m_dot_fg * (efg_comb - efg_exitvp);
StaExLosses  = (1 - eta_chimnex) * m_dot_fg*(efg_comb - ea);
ComExLosses  = (1 -  eta_combex) * m_dot_c * ec;

% Heat Exchanger (HE) losses and valve losses
if fh <= gamma
    HeatExLosses = - m_dot_vg*(state{ind1}.e - state{ind6}.e);
    for i = 1:fh
        indFR = i + 2;
        ind4i = alpha + beta*(i - 1) + 1; 
        HeatExLosses = HeatExLosses + FlowRates{2,indFR}*(state{ind4i}.e - state{ind8}.e);
    end
elseif fh > gamma
    HeatExLosses = - m_dot_vg*(state{ind1}.e - state{ind6BP}.e)...
                   - mv*(state{ind7BP}.e - state{ind6}.e);
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


%HeatExLosses = heatexchangerlosses(state,FlowRates)

%% Outputs
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
     'Heat Transfer Irreversibility in the Steam Generator'  'Heat Transfer Irreversibilities in FW heaters' ...
     'Stack Losses' 'Combustion Irreversibility'; P_el MecExLosses ConExLosses RotExLosses HeatExLosses ...
     TraExLosses StaExLosses ComExLosses};

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
function [ exergyFG ] = ExergyFlueGasses(tfg)
% Function Calculating the Exergy of the Flue Gasses temperature tfg in [°C]. 

% Old Version
% s0fg = EntropyFG(15);
% sfg  = EntropyFG(tfg);
% h0fg = EnthalpyFG(15);
% exergyFG = Exergy(h0fg,hf,s0fg,sfg)

delta_sfg  = EntropyFG(tfg) - EntropyFG(15);
delta_hfg  = cpfg(15,tfg)*(tfg - 15);

t0 = 273.15 + 15; %[K]
exergyFG = delta_hfg - t0*delta_sfg;
end

%% 
function [ er ] = ExergyR (ha_warm,t)
% Fonction calculating the exergy of the reagents with t the temperature of 
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
ea = cpair(15,Twarm)*(Twarm - 15) - t0*(EntropyAir(Twarm) - EntropyAir(15));
er = ea * lambda*ma1/(1 + lambda*ma1);
end

%%
function [ cpair ] = cpair (t1,t2)
% Function calculating the average specific heat for air between t1 and t2 in °C

T = linspace(273.15 + t1,273.15 + t2);
cpair = mean(0.21*janaf('c','O2',T) + 0.79*janaf('c','N2',T));
end

%%
function [ S ] = EntropyAir (ta)
% Function calculating the entropy of the flue gasses at temperature Tfg in °C


SO2  = janaf('s','O2',273.15 + ta);
SN2  = janaf('s','N2',273.15 + ta);
S    = (0.21*MWO2*SO2 + 0.79*MWN2*SN2)/MWair; %kJ / kgairK
end

%%
function [ cpfg ] = cpfg (T1,T2)
% Function calculating the average specific heat for flue gasses between 
% temperatures T1 and T2 in [°C]
   
T = linspace(273.15 + T1,273.15 + T2);
cpfg = mean(CO2_frac*janaf('c','CO2',T) + H2O_frac*janaf('c','H2O',T)...
         + O2_frac*janaf('c','O2',T)  +  N2_frac*janaf('c','N2',T));
end

%%
function [ tfg ] = TemperatureFG (hfg)
% Function calculating the adiabatic flame temperature of the flue gasses. 
% Take the combustion enthalpy as parameter in [kJ/kg]

tfg0 = 1500;                    % First guess for T combustion
cpfg0 = cpfg (0,tfg0);          % First guess for cp flue gasses

i = 0;
tfgprev = tfg0;
tfgnew  = hfg / cpfg0;
while (abs(tfgprev - tfgnew) > 1) && (i < 50)
    cpfgnew = cpfg (0,tfgnew);
    tfgprev = tfgnew;
    tfgnew  = hfg / cpfgnew;
    i = i + 1;
    if i >= 49
        warning ('max number of iterations in TemperatureFG reached')
    end
end
tfg  = tfgnew;               % OK has been verified 
    
end

%%
function [ S ] = EntropyFG (Tfg)
% Function calculating the entropy of the flue gasses at temperature Tfg in °C


SO2f  = janaf('s','O2',273.15 + Tfg);
SN2f  = janaf('s','N2',273.15 + Tfg);
SCO2f = janaf('s','CO2',273.15 + Tfg);
SH2Of = janaf('s','H2O',273.15 + Tfg);
S = CO2_frac*SCO2f + H2O_frac*SH2Of + O2_frac*SO2f +  N2_frac*SN2f; %kJ / kgfgK
end

%%
function flowrates = FRCreation(fr)
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

end

