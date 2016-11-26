function [ Output ] = SteamCycle( T_max, P_max, FH, RH )
%Steam Cycle
%   Input Arguments: 
%   - Max Steam Pressure [bar]
%   - Max Temperature [°C]
%   - Number of Feed Heaters [min value is 0]
%   - Reheating on or off
%   =======================================================================

%%
switch RH
    case 'off' 
        if FH == 0
            Output = RankineHirn(T_max,P_max);
        elseif FH > 0
            Output = FeedHeatersNoRH(T_max,P_max,FH);
        else
            warning('Negative number of feedheaters not allowed')
        end
    case 'on'
        if FH == 0
            Output = 2;
        elseif FH > 0
            Output = 3;
        else
            warning('Negative number of feedheaters not allowed')
            Output = 0;
        end
    otherwise
        warning('Unexpected reheating entry.')
        Output = 0;
end

%%
function [ Rankine ] = RankineHirn( t3, p3 )
% Rankine-Hirn: Simple Rankine-Hirn cycle.
%   Input Arguments: 
%   - Temperature at point 3 [°C]
%   - Pressure at point 3 [bar]
%   =======================================================================

n = 4;                     % number of states
state = State_creation(n,4,0); %Creation of the structure

% Point 3
state{3}.t = t3;
state{3}.p = p3;
state{3}.h = XSteam('h_pT',state{3}.p,state{3}.t);
state{3}.s = XSteam('s_pT',state{3}.p,state{3}.t);
state{3}.e = Exergy(state{3}.h,state{3}.s);

% Point 4 exit turbine
t4 = 32;
p4 = XSteam('psat_t',32);
eta_SiT = 0.88; % Isentropic Efficiency turbine
state{4} = turbine(state{3},t4,p4,eta_SiT);

% Point 1
state{1}.t = state{4}.t; %A rafiner peut-etre
state{1}.p = state{4}.p;
state{1}.x = 0;
state{1}.h = XSteam('hL_T',state{1}.t);
state{1}.s = XSteam('sL_t',state{1}.t);
state{1}.e = Exergy(state{1}.h,state{1}.s);

% Point 2 exit FWPump
eta_FWP  = 0.85; % Efficiency feedwater pump
kpdgen  = 0.90; % Pressure drop coefficient at steam generator. Determined via numerical examples in book
p2 = state{3}.p/kpdgen;
state{2} = pump(state{1},state{2},p2,eta_FWP);
Rankine = state;
end

%% 
function [ FHNoRH ] =  FeedHeatersNoRH(t3,p3,fh)
% Cycle with feedheaters but no reheating
%   Input Arguments: 
%   - Temperature at point 3 [°C]
%   - Pressure at point 3 [bar]
%   - Number of Feedheaters
%   =======================================================================

alpha = 8;                 %Number of fixed points of the schematic
beta  = 4;                 %Number of points per state: eg: beta = 3 means 4.1 6.1 7.1 or 4.2 6.2 7.2 
                           %It's in case we need to add one to take
                           %something more into account. eg before and
                           %after isenthalpic valves
gamma = 4;                 %Param for bache calculation
bache_pressure = 4.6;      %Param fixing the bache pressure

n = alpha + beta * fh;
state = State_creation(n,alpha,beta); %Creation of the structure

base = RankineHirn(t3,p3);
state{5} = base{1};
state{5}.States = '5';
state{3} = base{3};
state{4} = base{4};

% Les Points 4i
for i = 1:fh
    ind = alpha + beta*(i-1) + 1;
    state{ind}.h = state{4}.h + (state{3}.h - state{4}.h)*i / (fh + 1);
    state{ind}.s = state{3}.s + (state{4}.s - state{3}.s)*(fh + 1 - i) / (fh + 1);
    state{ind}.t = XSteam('T_hs',state{ind}.h,state{ind}.s);
    state{ind}.p = XSteam('p_hs',state{ind}.h,state{ind}.s);
    state{ind}.e = Exergy(state{ind}.h,state{ind}.s);
    if state{ind}.h < XSteam('hV_T',state{ind}.t)
        state{ind}.x = XSteam('x_ph',state{ind}.p,state{ind}.h);
    end
end

% Point 1 entry FWPump
state{1}.t = XSteam('Tsat_p',state{n-beta+1}.p) - 5;  %FWPump entry temp is sat temp of last bleed - 5°C.
state{1}.p = XSteam('psat_T',state{1}.t) + 10;   %FWPump entry pressure above saturation pressure.
state{1}.h = XSteam('h_pT',state{1}.p,state{1}.t);
state{1}.s = XSteam('s_pT',state{1}.p,state{1}.t);
state{1}.e = Exergy(state{1}.h,state{1}.s);

% Point 2 exit FWPump
eta_FWP  = 0.85; % Efficiency feedwater pump
kpdgen  = 0.90; % Pressure drop coefficient at steam generator. Determined via numerical examples in book
p2 = state{3}.p/kpdgen;
state{2} = pump(state{1},state{2},p2,eta_FWP);


position_bache = 1;               %References the roman number of the schematic
index_bleed_bache = alpha + 1;    %Index of the bleed entering the bache at higher pressure
if fh > gamma                     %Finding the bache if the number of bleeds is higher than gamma
    
    while state{index_bleed_bache}.p < 4.6
        position_bache = position_bache + 1;
        index_bleed_bache = index_bleed_bache + beta;
    end
    index_exit_bache = alpha + beta*position_bache; %References the exit of the bache 

    
    
    for i = 1:fh
        ind7i = alpha + beta*i -1;
        ind8i = alpha + beta*i;
        ind_bleed = alpha + beta*(i-1) + 1;
    
        if i == position_bache
            % Bache exit conditions point 8.bache = 7.bache because no valve
            state{index_exit_bache}.p = bache_pressure;          % Fixed pressure inside the bache
            state{index_exit_bache}.t = XSteam('Tsat_p',state{index_exit_bache}.p);
            state{index_exit_bache}.h = XSteam('hL_p',state{index_exit_bache}.p);
            state{index_exit_bache}.s = XSteam('sL_p',state{index_exit_bache}.p);
            state{index_exit_bache}.x = 0;
            state{index_exit_bache}.e = Exergy(state{index_exit_bache}.h,state{index_exit_bache}.s);
            
            str                              = state{index_exit_bache-1}.States;
            state{index_exit_bache-1}        = state{index_exit_bache};
            state{index_exit_bache-1}.States = str;
        else
            % Points 7.i exit of the bleed condensors before valves
            state{ind7i}.p = state{ind_bleed}.p;
            state{ind7i}.t = XSteam('Tsat_p',state{ind7i}.p);
            state{ind7i}.x = 0;
            state{ind7i}.h = XSteam('hL_p',state{ind7i}.p);
            state{ind7i}.s = XSteam('sL_p',state{ind7i}.p);
            state{ind7i}.e = Exergy(state{ind7i}.h,state{ind7i}.s);
             
             
            % Points 8.i exit of isenthalpic valves of the bleed condensors
            state{ind8i}.p = state{ind_bleed - beta}.p;
            state{ind8i}.t = state{ind7i}.t;
            state{ind8i}.h = state{ind7i}.h;
            state{ind8i}.x = XSteam('x_ph',state{ind8i}.p,state{ind8i}.h);
            state{ind8i}.s = XSteam('s_ph',state{ind8i}.p,state{ind8i}.h);
            state{ind8i}.e = Exergy(state{ind8i}.h,state{ind8i}.s);
        end
    end
    
    % Pressure point 6.0 (exit Condensor Pump)
    state{6}.p = state{index_exit_bache}.p;
    
else  %Without bache
    state{6}.p = state{1}.p;
    
    % Points 8.i exit of isenthalpic valves of the bleed condensers
    for i = 1:fh
        ind7i = alpha + beta*i -1;
        ind8i = alpha + beta*i;
        ind_bleed = alpha + beta*(i-1) + 1;

        % Points 7.i exit of the bleed condensors before valves
        state{ind7i}.p = state{ind_bleed}.p;
        state{ind7i}.t = XSteam('Tsat_p',state{ind7i}.p);
        state{ind7i}.x = 0;
        state{ind7i}.h = XSteam('hL_p',state{ind7i}.p);
        state{ind7i}.s = XSteam('sL_p',state{ind7i}.p);
        state{ind7i}.e = Exergy(state{ind7i}.h,state{ind7i}.s);


        % Points 8.i exit of isenthalpic valves of the bleed condensors
        state{ind8i}.p = state{ind_bleed - beta}.p;
        state{ind8i}.t = state{ind7i}.t;
        state{ind8i}.h = state{ind7i}.h;
        state{ind8i}.x = XSteam('x_ph',state{ind8i}.p,state{ind8i}.h);
        state{ind8i}.s = XSteam('s_ph',state{ind8i}.p,state{ind8i}.h);
        state{ind8i}.e = Exergy(state{ind8i}.h,state{ind8i}.s);
    end
end

eta_CP = 0.85;
state{6} = pump(state{5},state{6},state{6}.p,eta_CP);

% Point 7 exit Sub-Cooler
state{7}.t = state{6}.t + 5;
state{7}.p = state{(alpha + beta - 1)}.p;
state{7}.h = XSteam('h_pT',state{7}.p,state{7}.t);
state{7}.s = XSteam('s_pT',state{7}.p,state{7}.t);
state{7}.e = Exergy(state{7}.h,state{7}.s);

% Point 8 global bleedwater at entrance condensor
state{8}.p = state{4}.p;
state{8}.t = state{7}.t;
state{8}.h = state{7}.h;
state{8}.x = XSteam('x_ph',state{8}.p,state{8}.h);
state{8}.s = XSteam('s_ph',state{8}.p,state{8}.h);
state{8}.e = Exergy(state{8}.h,state{8}.s);



FHNoRH = state;
end

%% 
function ex = Exergy (h,s)
    t0 = 273.15 + 15; %[K]
    p0 = 0.01704;
    h0 = 63.0;
    s0 = 0.224;
    ex = h - h0 - t0*(s - s0);
end

function muT = muT(t)
% Function to extrapolate the value of muT of water at different temperatures
% using values of the LMECA1855 exercice book at 30 bar. Pressure has not a
% big effect on muT
% Input Variables:
%   - Temperature of the water. Max temp = 230 °C
% =========================================================================
   
    T   = [0.02 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230];
    MuT = [0.10140 0.09720 0.09385 0.09101 0.08847 0.08614 0.08392 0.08175 ...
           0.07960 0.07743 0.07520 0.07288 0.07044 0.06784 0.06504 0.06300 ...
           0.05865 0.05494 0.05079 0.04610 0.04073 0.03454 0.02730 0.01871];
    
    ind = 1;
    while t > T(ind) %ind will be the index of the closest higher temp of the test temp. 
        ind = ind + 1;
    end
muT = (MuT(ind) - MuT(ind-1)) * (t - T(ind-1)) / 10 + MuT(ind-1);

end

function Output = pump(state_in,state_out,p,eta)
% Function calculating the output state of a pump.
% Input Variables:
%   - State at entrance of the pump
%   - Desired exit pressure
%   - Efficiency of the pump
% =========================================================================
    
v_LH2O   = 0.001005; %[m³/kg] volume massique de l'eau

state_out.States = state_out.States;
state_out.p = p;
state_out.h = v_LH2O * (state_out.p - state_in.p)*100 / eta + state_in.h; % *e2 pour obtenir kJ/kg (si e5 on obtient des joules..)
muT            = muT(state_in.t);                                           % muT approx with temp before pump
cp2            = (XSteam('Cp_pT',state_out.p,state_in.t) + XSteam('CpL_p',state_in.p)) / 2; 
state_out.t = state_in.t + (v_LH2O*100 - muT)*(state_out.p - state_in.p)/ (cp2);
state_out.s = XSteam ('s_pT',state_out.p,state_out.t);
state_out.e = Exergy(state_out.h,state_out.s);

Output = state_out;
end

function Output = turbine (state,t_low,p_low,eta)
% Function calculating the output state of a turbine.
% Input Variables:
%   - State at entrance of the turbine
%   - Exit temperature
%   - Exit pressure
%   - Isentropic efficiency of the turbine
% =========================================================================
   
    state_out = State_creation(1,1,1);
    
    state_out{1}.States = state.States + 1;
    state_out{1}.t = t_low;
    state_out{1}.p = p_low;
    s_isos  = state.s;
    
    t_sat_test = XSteam('Tsat_p',p_low);
    if t_low < t_sat_test
        x_isos  = XSteam('x_ps',p_low,s_isos); % - XSteam('sL_T',state_out{1}.t)) / (XSteam('sV_T',state_out{1}.t) - XSteam('sL_T',state_out{1}.t));
        h_isos  = XSteam('h_px',p_low,x_isos); % x4s*XSteam('hV_T',state_out{1}.t) + (1-x4s)*XSteam('hL_T',state_out{1}.t);
    else
        h_isos  = XSteam('h_ps',p_low,s_isos);
    end

    state_out{1}.h = eta*(h_isos - state.h) + state.h;
    state_out{1}.x = XSteam('x_ph',state_out{1}.p,state_out{1}.h);
    state_out{1}.s = state_out{1}.x*XSteam('sV_T',state_out{1}.t) + (1 - state_out{1}.x)*XSteam('sL_T',state_out{1}.t);
    state_out{1}.e = Exergy(state_out{1}.h,state_out{1}.s);
    
    Output =  state_out{1};
end

end
