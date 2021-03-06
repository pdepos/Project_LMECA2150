function [ Output ] = SCStates( T_max, T_min, P_max, FH, RH )
%Steam Cycle
%   Input Arguments: 
%   - Max Steam Pressure [bar]
%   - Max Temperature [�C]
%   - Number of Feed Heaters [min value is 0]
%   - Reheating on or off
%   =======================================================================

%%
step = 10000;
switch RH
    case 'off' 
        if (T_max > 700) || (T_max < 200) || (P_max > 100) || (P_max < 10)
            warning('The cycle you are trying to test is not reasonable. Please verify pressure or temperature inputs')
            return
        end
        r = 0;
        if FH == 0
            Output = RankineHirn(T_max, T_min, P_max);
        elseif FH > 0
            Output = BasicFeedHeating(T_max, T_min, P_max, FH, r);
        else
            warning('Negative number of feedheaters not allowed')
        end
    case 'on'
        r = 1;
        if (T_max > 700) || (T_max < 200) || (P_max > 400) || (P_max < 50)
            warning('The cycle you are trying to test is not reasonable. Please verify pressure or temperature inputs')
            return
        end
        if FH == 0
            Output = BasicReHeating(T_max, T_min, P_max);
        elseif FH > 0
            Output = ReHeatingAndFH(T_max, T_min, P_max, FH, r);
        else
            warning('Negative number of feedheaters not allowed')
            Output = 0;
        end
    otherwise
        warning('Unexpected reheating entry.')
        Output = 0;
end

%%
function [ Rankine ] = RankineHirn( t_max, t_min, p_max )
% Rankine-Hirn: Simple Rankine-Hirn cycle.
%   Input Arguments: 
%   - Temperature at point 3 [�C]
%   - Pressure at point 3 [bar]
%   =======================================================================

n = 6;                         % number of states
state = State_creation(n,6,0); % Creation of the structure

ind2i  = 3;
ind2ii = 4;
ind3   = 5;
ind4   = 6;
eta_FWP = 0.85; % Efficiency feedwater pump
eta_SiT = 0.88; % Isentropic Efficiency turbine
kpdgen  = 1.10; % Pressure drop coefficient at steam generator. Determined via numerical examples in book

% Point 3
state{ind3}.t = t_max;
state{ind3}.p = p_max;
state{ind3}.h = XSteam('h_pT',state{ind3}.p,state{ind3}.t);
state{ind3}.s = XSteam('s_pT',state{ind3}.p,state{ind3}.t);
state{ind3}.e = Exergy(state{ind3}.h,state{ind3}.s);

% Point 4 exit turbine
state{ind4}.p = XSteam('psat_T',t_min);
state{ind4} = turbine(state{ind3},state{ind4},eta_SiT);

% Point 1
state{1}.t = state{ind4}.t - 3; %A rafiner peut-etre
state{1}.p = XSteam('Psat_T',state{1}.t); %state{ind4}.p;
state{1}.x = 0;
state{1}.h = XSteam('hL_T',state{1}.t);
state{1}.s = XSteam('sL_t',state{1}.t);
state{1}.e = Exergy(state{1}.h,state{1}.s);

% Point 2 exit FWPump
state{2}.p = state{ind3}.p*kpdgen;                 
state{2} = pump(state{1},state{2},eta_FWP);

% Point 2'
state{ind2i}.p = state{ind3}.p;
state{ind2i}.t = XSteam('Tsat_p',state{ind2i}.p);
state{ind2i}.x = 0;
state{ind2i}.h = XSteam('hL_p',state{ind2i}.p);
state{ind2i}.s = XSteam('sL_p',state{ind2i}.p);
state{ind2i}.e = Exergy(state{ind2i}.h,state{ind2i}.s);

% Point 2''
state{ind2ii}.p = state{ind3}.p;
state{ind2ii}.t = XSteam('Tsat_p',state{ind2ii}.p);
state{ind2ii}.x = 1;
state{ind2ii}.h = XSteam('hV_p',state{ind2ii}.p);
state{ind2ii}.s = XSteam('sV_p',state{ind2ii}.p);
state{ind2ii}.e = Exergy(state{ind2ii}.h,state{ind2ii}.s);

Rankine = state;
end

%% 
function [ BFH ] =  BasicFeedHeating( t_max, t_min, p_max, fh , rh_param)
% Cycle with feedheaters but no reheating
%   Input Arguments: 
%   - Temperature at point 3 [�C]
%   - Pressure at point 3 [bar]
%   - Number of Feedheaters
%   =======================================================================

alpha = 10;                %Number of fixed points of the schematic
beta  = 4;                 %Number of points per state: eg: beta = 3 means 4.1 6.1 7.1 or 4.2 6.2 7.2 
                           %It's in case we need to add one to take
                           %something more into account. eg before and
                           %after isenthalpic valves
gamma    = 4;              %Param for bache calculation
p_bache0 = 2;            %Param fixing the bache pressure
eta_FWP  = 0.85;           %Efficiency feedwater pump
eta_BP   = 0.85;           %Efficiency of the Bache Pump
eta_CP   = 0.85;           %Efficiency Condenser Pump
eta_SiT  = 0.88;
kpdgen   = 1.10;           %Pressure drop coefficient at steam generator. 
                           %Determined via numerical examples in book
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

n = alpha + beta * fh;
state = State_creation(n,alpha,beta); %Creation of the structure

base = RankineHirn(t_max,t_min,p_max);
state{ind2i}  = base{ind2i};
state{ind2ii} = base{ind2ii};
state{ind3}   = base{ind3};
state{ind4}   = base{ind4};
state{ind5}   = base{ind1};
state{ind5}.States = '5';


% Les Points 4i
    p34   = linspace(state{ind3}.p,state{ind4}.p,step);
    s34_s = linspace(state{ind3}.s,state{ind3}.s,step);
    h34_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34,s34_s);
    h34   = zeros(1,step);
    for i = 1:step
        h34(i) = state{ind3}.h + eta_SiT*(h34_s(i) - state{ind3}.h);
    end
    s34   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34,h34);
    
for i = 1:fh
    ind4i = alpha + beta*(i-1) + 1;
    state{ind4i}.h = state{ind4}.h + (state{ind3}.h - state{ind4}.h)*i / (fh + 1);
    error = abs(state{ind4i}.h - h34(1));
    j = 1;
    while (error >= 2) && (j < step)
        j = j + 1;
        error = abs(state{ind4i}.h - h34(j));
        if j == 9999
            warning('j = 9999')
        end
    end
    state{ind4i}.s = s34(j);
    state{ind4i}.t = XSteam('T_hs',state{ind4i}.h,state{ind4i}.s);
    state{ind4i}.p = XSteam('p_hs',state{ind4i}.h,state{ind4i}.s);
    state{ind4i}.e = Exergy(state{ind4i}.h,state{ind4i}.s);
    if state{ind4i}.h < XSteam('hV_T',state{ind4i}.t)
        state{ind4i}.x = XSteam('x_ph',state{ind4i}.p,state{ind4i}.h);
    end
end

% Point 1 entry FWPump
if rh_param == 1;                            %If used for RH Cycle, FWP entry Temp depends on exit pressure HP Turbine
    p_4HP = p_max * kpdgen;                  %Pressure exit HP Turbine if RH Cycle
    state{ind1}.t = XSteam('Tsat_p',p_4HP);  %To consider the desuperheaters in a simple way we decided to fix 
                                             %the temperature of the FW at the entry of the FWPump to be equal to 
                                             %the sat temp of the HP Turbine Bleed in approximation to what we saw in the book.
    state{ind1}.p = XSteam('psat_T',state{ind1}.t) + 10;        %FWPump entry pressure above saturation pressure.
    state{ind1}.h = XSteam('h_pT',state{ind1}.p,state{ind1}.t);
    state{ind1}.s = XSteam('s_pT',state{ind1}.p,state{ind1}.t);
    state{ind1}.e = Exergy(state{ind1}.h,state{ind1}.s);
else
    state{ind1}.t = XSteam('Tsat_p',state{n-beta+1}.p) - 5;     %FWPump entry temp is sat temp of last bleed - 5�C.
    state{ind1}.p = XSteam('psat_T',state{ind1}.t) + 2;         %FWPump entry pressure above saturation pressure.
    state{ind1}.h = XSteam('h_pT',state{ind1}.p,state{ind1}.t);
    state{ind1}.s = XSteam('s_pT',state{ind1}.p,state{ind1}.t);
    state{ind1}.e = Exergy(state{ind1}.h,state{ind1}.s);
end

% Point 2 exit FWPump
state{ind2}.p = state{ind3}.p * kpdgen;
state{ind2} = pump(state{ind1},state{ind2},eta_FWP);


position_bache = 1;               %References the roman number of the schematic
index_bleed_bache = alpha + 1;    %Index of the bleed entering the bache at higher pressure
if fh > gamma                     %Finding the bache if the number of bleeds is higher than gamma
    while (state{index_bleed_bache}.p < p_bache0) && (position_bache <= fh)
        position_bache = position_bache + 1;
        index_bleed_bache = index_bleed_bache + beta;
    end
    index_exit_bache = alpha + beta*position_bache; %References the exit 8bache of the bache 
    p_bache = state{index_bleed_bache}.p;
    
    for i = 1:fh
        ind7i = alpha + beta*i - 1;
        ind8i = alpha + beta*i;
        ind_bleed = alpha + beta*(i-1) + 1;
        
        if i == position_bache
            % Bache exit conditions point 8.bache = 7.bache because no valve
            state{index_exit_bache}.p = p_bache;          % Fixed pressure inside the bache
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
            if i == 1
                % Points 8.1 same state as 7.1
                state{ind8i}.p = state{ind7i}.p;
                state{ind8i}.t = state{ind7i}.t;
                state{ind8i}.h = state{ind7i}.h;
                state{ind8i}.x = XSteam('x_ph',state{ind8i}.p,state{ind8i}.h);
                state{ind8i}.s = XSteam('s_ph',state{ind8i}.p,state{ind8i}.h);
                state{ind8i}.e = Exergy(state{ind8i}.h,state{ind8i}.s);
            else 
                % Points 8.i exit of isenthalpic valves
                state{ind8i}.p = state{ind_bleed - beta}.p;
                state{ind8i}.t = state{ind7i}.t;
                state{ind8i}.h = state{ind7i}.h;
                state{ind8i}.x = XSteam('x_ph',state{ind8i}.p,state{ind8i}.h);
                state{ind8i}.s = XSteam('s_ph',state{ind8i}.p,state{ind8i}.h);
                state{ind8i}.e = Exergy(state{ind8i}.h,state{ind8i}.s);
            end
        end        
    end
    
    % Point 6.0 (exit Condensor Pump)
    state{ind6}.p = state{index_exit_bache}.p;
    state{ind6} = pump(state{ind5},state{ind6},eta_CP);
    
    % Points 6.i FW circuit except the first one at this stage
    for i = 2:fh
        ind6i = alpha + beta*i - 2;
        ind7i = alpha + beta*(i - 1) - 1;  % point 7.(i-1)
        
        if i <= position_bache              % Points 6.i before bache
            state{ind6i}.p = state{ind6}.p;
            state{ind6i}.t = state{ind7i}.t - 5;
            state{ind6i}.h = XSteam('h_pT',state{ind6i}.p,state{ind6i}.t);
            state{ind6i}.s = XSteam('s_pT',state{ind6i}.p,state{ind6i}.t);
            state{ind6i}.e = Exergy(state{ind6i}.h,state{ind6i}.s);
        elseif (i == position_bache + 1)   % Point 6 behind Bache Pump
            state{ind6i}.p = state{ind1}.p;
            state{ind6i} = pump(state{ind7i},state{ind6i},eta_BP);
        else                               % Points 6.i behind the bache pump
            state{ind6i}.p = state{ind1}.p;
            state{ind6i}.t = state{ind7i}.t - 5;
            state{ind6i}.h = XSteam('h_pT',state{ind6i}.p,state{ind6i}.t);
            state{ind6i}.s = XSteam('s_pT',state{ind6i}.p,state{ind6i}.t);
            state{ind6i}.e = Exergy(state{ind6i}.h,state{ind6i}.s);
        end
    end
    
else  %Without bache
    for i = 1:fh
        ind7i = alpha + beta*i - 1;
        ind8i = alpha + beta*i;
        ind_bleed = alpha + beta*(i-1) + 1;

        % Points 7.i exit of the bleed condensors before valves
        state{ind7i}.p = state{ind_bleed}.p;
        state{ind7i}.t = XSteam('Tsat_p',state{ind7i}.p);
        state{ind7i}.x = 0;
        state{ind7i}.h = XSteam('hL_p',state{ind7i}.p);
        state{ind7i}.s = XSteam('sL_p',state{ind7i}.p);
        state{ind7i}.e = Exergy(state{ind7i}.h,state{ind7i}.s);
        if i == 1
            % Points 8.1 same state as 7.1
            state{ind8i}.p = state{ind7i}.p;
            state{ind8i}.t = state{ind7i}.t;
            state{ind8i}.h = state{ind7i}.h;
            state{ind8i}.x = XSteam('x_ph',state{ind8i}.p,state{ind8i}.h);
            state{ind8i}.s = XSteam('s_ph',state{ind8i}.p,state{ind8i}.h);
            state{ind8i}.e = Exergy(state{ind8i}.h,state{ind8i}.s);
        else 
            % Points 8.i exit of isenthalpic valves
            state{ind8i}.p = state{ind_bleed - beta}.p;
            state{ind8i}.t = state{ind7i}.t;
            state{ind8i}.h = state{ind7i}.h;
            state{ind8i}.x = XSteam('x_ph',state{ind8i}.p,state{ind8i}.h);
            state{ind8i}.s = XSteam('s_ph',state{ind8i}.p,state{ind8i}.h);
            state{ind8i}.e = Exergy(state{ind8i}.h,state{ind8i}.s);
        end
    end
            
    % Point 6.0 (exit Condensor Pump) 
    state{ind6}.p = state{ind1}.p;
    state{ind6} = pump(state{ind5},state{ind6},eta_CP);
    
    % Points 6.i FW circuit except the first one at this stage
    for i = 2:fh
        ind6i = alpha + beta*i - 2;
        ind7i = alpha + beta*(i - 1) - 1;  % point 7.(i-1)
        
        state{ind6i}.p = state{ind6}.p;
        state{ind6i}.t = state{ind7i}.t - 5;
        state{ind6i}.h = XSteam('h_pT',state{ind6i}.p,state{ind6i}.t);
        state{ind6i}.s = XSteam('s_pT',state{ind6i}.p,state{ind6i}.t);
        state{ind6i}.e = Exergy(state{ind6i}.h,state{ind6i}.s);
    end
    
end

% Point 7 exit Sub-Cooler
state{ind7}.t = state{ind6}.t + 5;
state{ind7}.p = state{(alpha + beta - 1)}.p;
state{ind7}.h = XSteam('h_pT',state{ind7}.p,state{ind7}.t);
state{ind7}.s = XSteam('s_pT',state{ind7}.p,state{ind7}.t);
state{ind7}.e = Exergy(state{ind7}.h,state{ind7}.s);

% Point 8 global bleedwater at entrance condensor
state{ind8}.p = state{ind4}.p;
state{ind8}.t = state{ind7}.t;
state{ind8}.h = state{ind7}.h;
state{ind8}.x = XSteam('x_ph',state{ind8}.p,state{ind8}.h);
state{ind8}.s = XSteam('s_ph',state{ind8}.p,state{ind8}.h);
state{ind8}.e = Exergy(state{ind8}.h,state{ind8}.s);


BFH = state;
end

%%
function [ RHNoFH ] = BasicReHeating( t_max, t_min, hp )
  
ind1    = 1;
ind2    = 2;
ind2i   = 3;
ind2ii  = 4;
ind3HP  = 5;
ind4HP  = 6;
ind3LP  = 7;
ind4LP  = 8;

p_3LP     = 0.12*hp;
kpdgen    = 1.10;
eta_SiTHP = 0.90; % Isentropic Efficiency HP Turbine
eta_SiTLP = 0.88; % Isentropic Efficiency LP Turbine
eta_FWP   = 0.85;

n = 8;
state = State_creation(n,n,0);

% Point 3HP
state{ind3HP}.p = hp;
state{ind3HP}.t = t_max;
state{ind3HP}.h = XSteam('h_pT',state{ind3HP}.p,state{ind3HP}.t);
state{ind3HP}.s = XSteam('s_pT',state{ind3HP}.p,state{ind3HP}.t);
state{ind3HP}.e = Exergy(state{ind3HP}.h,state{ind3HP}.s);

% Point 4HP
state{ind4HP}.p = p_3LP * kpdgen;
state{ind4HP}   = turbine(state{ind3HP},state{ind4HP},eta_SiTHP);

% Point 3LP
state{ind3LP}.p = p_3LP;
state{ind3LP}.t = t_max;
state{ind3LP}.h = XSteam('h_pT',state{ind3LP}.p,state{ind3LP}.t);
state{ind3LP}.s = XSteam('s_pT',state{ind3LP}.p,state{ind3LP}.t);
state{ind3LP}.e = Exergy(state{ind3LP}.h,state{ind3LP}.s);

% Point 4LP
state{ind4LP}.p = XSteam('Psat_T',t_min);
state{ind4LP} = turbine(state{ind3LP},state{ind4LP},eta_SiTLP);

% Point 1
state{ind1}.t = state{ind4LP}.t - 3;
state{ind1}.p = XSteam('Psat_T',state{1}.t);
state{ind1}.x = 0;
state{ind1}.h = XSteam('hL_T',state{ind1}.t);
state{ind1}.s = XSteam('sL_t',state{ind1}.t);
state{ind1}.e = Exergy(state{ind1}.h,state{ind1}.s);

% Point 2 exit FWPump
state{ind2}.p = state{ind3HP}.p * kpdgen;
state{ind2} = pump(state{ind1},state{ind2},eta_FWP);

% Point 2'
state{ind2i}.p = state{ind3HP}.p;
state{ind2i}.t = XSteam('Tsat_p',state{ind2i}.p);
state{ind2i}.x = 0;
state{ind2i}.h = XSteam('hL_p',state{ind2i}.p);
state{ind2i}.s = XSteam('sL_p',state{ind2i}.p);
state{ind2i}.e = Exergy(state{ind2i}.h,state{ind2i}.s);

% Point 2''
state{ind2ii}.p = state{ind3HP}.p;
state{ind2ii}.t = XSteam('Tsat_p',state{ind2ii}.p);
state{ind2ii}.x = 1;
state{ind2ii}.h = XSteam('hV_p',state{ind2ii}.p);
state{ind2ii}.s = XSteam('sV_p',state{ind2ii}.p);
state{ind2ii}.e = Exergy(state{ind2ii}.h,state{ind2ii}.s);

RHNoFH = state;
end

%%
function [ RHandFH ] = ReHeatingAndFH( t_max, t_min, hp, fh, rh_param )
alpha = 12;                %Number of fixed points of the schematic
beta  = 4;                 %Number of points per state: eg: beta = 3 means 4.1 6.1 7.1 or 4.2 6.2 7.2 
                           %It's in case we need to add one to take
                           %something more into account. eg before and
                           %after isenthalpic valves
alpha_FH_base = 10;        %Fixed points in simple FH function to get the indexes in the FH_base
eta_FWP   = 0.85;
kpdgen    = 1.10;          %Pressure drop coefficient at steam generator. 
                           %Determined via numerical examples in book
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

n = alpha + beta * fh + 4;            % + 4 because we consider the HP bleed outside of the LP bleeds. 
state = State_creation(n,alpha,beta); %Creation of the structure

% Base Points of the ReHeating Cycle
RH_base = BasicReHeating(t_max,t_min,hp);
state{ind2i}  = RH_base{ind2i};
state{ind2ii} = RH_base{ind2ii};
state{ind3HP} = RH_base{ind3HP};
state{ind4HP} = RH_base{ind4HP};
state{ind3LP} = RH_base{ind3LP};
state{ind4LP} = RH_base{ind4LP};

% Points of the Feed Heating Cycle
FH_base = BasicFeedHeating(t_max,t_min,state{ind3LP}.p,fh,rh_param);
state{ind1} = FH_base{ind1};
state{ind2}.p = hp * kpdgen;
state{ind2} = pump(state{ind1},state{ind2},eta_FWP);

state{ind5} = FH_base{ind5-2};
state{ind6} = FH_base{ind6-2};
state{ind7} = FH_base{ind7-2};
state{ind8} = FH_base{ind8-2};

for i = 1:fh
    % Points 4.i
    ind4i = alpha + beta*(i-1) + 1;
    ind4i_base = alpha_FH_base + beta*(i-1) + 1;
    state{ind4i} = FH_base{ind4i_base};
    
    % Points 6.i
    ind6i = alpha + beta*(i-1) + 2;
    ind6i_base = alpha_FH_base + beta*(i-1) + 2;
    state{ind6i} = FH_base{ind6i_base};
    
    % Points 7.i 
    ind7i = alpha + beta*(i-1) + 3;
    ind7i_base = alpha_FH_base + beta*(i-1) + 3;
    state{ind7i} = FH_base{ind7i_base};
    
    % Points 8.i
    ind8i = alpha + beta*(i-1) + 4;
    ind8i_base = alpha_FH_base + beta*(i-1) + 4;
    state{ind8i} = FH_base{ind8i_base};
end

% Points around the HP Bleed
% Point 4HP 
str1 = state{n - beta + 1}.States;
str2 = '=';
str3 = state{ind4HP}.States;
state{n - beta + 1} = state{ind4HP};
state{n - beta + 1}.States = [str1 str2 str3];

% Point 7HP
state{n - beta + 3}.p = state{n - beta + 1}.p;
state{n - beta + 3}.t = XSteam('Tsat_p',state{n - beta + 3}.p);
state{n - beta + 3}.x = 0;
state{n - beta + 3}.h = XSteam('hL_p',state{n - beta + 3}.p);
state{n - beta + 3}.s = XSteam('sL_p',state{n - beta + 3}.p);
state{n - beta + 3}.e = Exergy(state{n - beta + 3}.h,state{n - beta + 3}.s);

% Point 8HP
state{n - beta + 4}.p = state{n - 2*beta + 1}.p;
state{n - beta + 4}.t = state{n - beta + 3}.t;
state{n - beta + 4}.h = state{n - beta + 3}.h;
state{n - beta + 4}.x = XSteam('x_ph',state{n - beta + 4}.p,state{n - beta + 4}.h);
state{n - beta + 4}.s = XSteam('s_ph',state{n - beta + 4}.p,state{n - beta + 4}.h);
state{n - beta + 4}.e = Exergy(state{n - beta + 4}.h,state{n - beta + 4}.s);
  
% Point 6HP
state{n - beta + 2}.p = state{ind1}.p;
state{n - beta + 2}.t = state{n - 2*beta + 3}.t - 5;
state{n - beta + 2}.h = XSteam('h_pT',state{n - beta + 2}.p,state{n - beta + 2}.t);
state{n - beta + 2}.s = XSteam('s_pT',state{n - beta + 2}.p,state{n - beta + 2}.t);
state{n - beta + 2}.e = Exergy(state{n - beta + 2}.h,state{n - beta + 2}.s);


RHandFH = state;
end

%%
function ex = Exergy (h,s)
    t0 = 273.15 + 15; %[K]
    h0 = 62.96;
    s0 = 0.2244;
    ex = h - h0 - t0*(s - s0);
end

function Output = pump(state_in,state_out,eta)
% Function calculating the output state of a pump.
% Input Variables:
%   - State at entrance of the pump
%   - Desired exit pressure
%   - Efficiency of the pump
% =========================================================================
    
    v_LH2O   = 0.001005; %[m�/kg] volume massique de l'eau

    state_out.h = v_LH2O * (state_out.p - state_in.p)*100 / eta + state_in.h; % *e2 pour obtenir kJ/kg (si e5 on obtient des joules..)
    state_out.t = XSteam('T_ph',state_out.p,state_out.h);
    state_out.s = XSteam ('s_ph',state_out.p,state_out.h);
    state_out.e = Exergy(state_out.h,state_out.s);

    Output = state_out;
end

function Output = turbine (state_in,state_out,eta)
% Function calculating the output state of a turbine.
% Input Variables:
%   - State at entrance of the turbine. Defines entry entropy and enthalpy
%   - State at exit turbine. Defines exit pressure
%   - Isentropic efficiency of the turbine
% =========================================================================

    h_isos  = XSteam('h_ps',state_out.p,state_in.s);    
    state_out.h = eta*(h_isos - state_in.h) + state_in.h;   
    state_out.s = XSteam('s_ph',state_out.p,state_out.h);   
    state_out.t = XSteam('T_ph',state_out.p,state_out.h);
    state_out.e = Exergy(state_out.h,state_out.s);

    if state_out.t == XSteam('Tsat_p',state_out.p);
        state_out.x = XSteam('x_ph',state_out.p,state_out.h);
    else
        state_out.x = NaN;
    end

    Output =  state_out;
end

end
