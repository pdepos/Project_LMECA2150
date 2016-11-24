function [ Output ] = SteamCycle( T_max, P_max, FH, RH )
%Steam Cycle
%   Input Arguments: 
%   - Max Steam Pressure [bar]
%   - Max Temperature [�C]
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
%   - Temperature at point 3 [�C]
%   - Pressure at point 3 [bar]
%   =======================================================================

n = 4;                     % number of states
state = State_creation(n); %Creation of the structure

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
state{2} = pump(state{1},p2,eta_FWP);

Rankine = state;
end

%% 
function [ FHNoRH ] =  FeedHeatersNoRH(t3,p3,fh)
% Cycle with feedheaters but no reheating
%   Input Arguments: 
%   - Temperature at point 3 [�C]
%   - Pressure at point 3 [bar]
%   - Number of Feedheaters
%   =======================================================================

%   Note pour moi-m�me: Premi�rve version pour un nombre de FH bas, en ne consid�rant pas de
%   bache. Si FH > 4 alors on peut placer une bache, probablement entre le
%   3eme et 4eme FH.

n = 8 + 3 * fh;
state = State_creation(n); %Creation of the structure

base = RankineHirn(t3,p3);
state{5} = base{1};
%state{2} = base{2}; %ATTENTION cet etat n'est pas bon
state{3} = base{3};
state{4} = base{4};

% Les Points 4i
for i = 1:fh
    ind = 8 + 3*(i-1) + 1;
    state{ind}.h = state{4}.h + (state{3}.h - state{4}.h)*i / (fh + 1);
    state{ind}.s = state{3}.s + (state{4}.s - state{3}.s)*(fh + 1 - i) / (fh + 1);
    state{ind}.t = XSteam('T_hs',state{ind}.h,state{ind}.s);
    state{ind}.p = XSteam('p_hs',state{ind}.h,state{ind}.s);
    state{ind}.e = Exergy(state{ind}.h,state{ind}.s);
    if state{ind}.h < XSteam('hV_T',state{ind}.t)
        state{ind}.x = XSteam('x_ph',state{ind}.p,state{ind}.h);
                   % (state{ind}.h - XSteam('hL_T',state{ind}.t)) / ...
                   % (XSteam('hV_T',state{ind}.t) - XSteam('hL_T',state{ind}.t));
    end
end

% Point 1, entry FWPump
state{1}.t = XSteam('Tsat_p',state{n-2}.p) - 5;  %FWPump entry temp is sat temp of last bleed - 5�C.
state{1}.p = XSteam('psat_T',state{1}.t) + 10;   %FWPump entry pressure above saturation pressure.
state{1}.h = XSteam('h_pT',state{1}.p,state{1}.t);
state{1}.s = XSteam('s_pT',state{1}.p,state{1}.t);
state{1}.e = Exergy(state{1}.h,state{1}.s);

% Point 2, exit FWPump
eta_FWP  = 0.85; % Efficiency feedwater pump
kpdgen  = 0.90; % Pressure drop coefficient at steam generator. Determined via numerical examples in book
p2 = state{3}.p/kpdgen;
state{2} = pump(state{1},p2,eta_FWP);


position_before_bache = 1;
index_before_bache = 8 + 1;
while state{index_before_bache}.p < 2
    position_before_bache_new = position_before_bache + 1;
    position_before_bache = position_before_bache_new;
    index_before_bache_new = 8 + 3*(position_before_bache - 1) + 1;
    index_before_bache = index_before_bache_new;
end
position_bache = position_before_bache + 1


% Point 6.0 
% Note pour Rapport: On fixe la pression � mi-chemin entre la pression de
% sortie des deux feedheaters qui entourent la bache. Dans notre cycle on a
% toujours une bache pour avoir de la vapeur de qualit�. Elle se trouve 
% parfois en premi�re position lorsque le nombre de soutirages est petit.

%state{6}.p = (state{8 + 3*(position_bache - 1) + 1}.p + state{8 + 3*(position_before_bache - 1) + 1}.p)/2;






% Conversion de Structure vers matrice pour l'affichage
FHNoRH = ones(n:7);
for i = 1:n
        FHNoRH(i,:)    = [i state{i}.t state{i}.p state{i}.x ...
                              state{i}.h state{i}.s state{i}.e];
end
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

function Output = pump(state,p,eta)
    v_LH2O   = 0.001005; %[m�/kg] volume massique de l'eau
    state_out = State_creation(1);
    
    state_out{1}.p = p;
    state_out{1}.h = v_LH2O * (state_out{1}.p - state.p)*100 / eta + state.h; % *e2 pour obtenir kJ/kg (si e5 on obtient des joules..)
    muT            = muT(state.t);                                           % muT approx with temp before pump
    cp2            = (XSteam('Cp_pT',state_out{1}.p,state.t) + XSteam('CpL_p',state.p)) / 2; 
    state_out{1}.t = state.t + (v_LH2O*100 - muT)*(state_out{1}.p - state.p)/ (cp2);
    state_out{1}.s = XSteam ('s_pT',state_out{1}.p,state_out{1}.t);
    state_out{1}.e = Exergy(state_out{1}.h,state_out{1}.s);
    
    Output = state_out{1};
end

function Output = turbine (state,t_low,p_low,eta)
    state_out = State_creation(1);
    
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
