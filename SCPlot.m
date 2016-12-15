function [] = SCPlot( state, RH, FH )
% Function plotting the TS & HS diagrams of the Steam Cycle.

ind1    = 1;
ind2    = 2;
ind2i   = 3;
ind2ii  = 4;
step2 = 500;
eta_SiTHP = 0.90;
eta_SiTLP = 0.88;

switch RH
    case 'off' 
        if FH == 0
           alpha = 6;
           ind3  = 5;
           ind4  = 6;
        elseif FH > 0
           alpha = 10;
           beta  = 4;
           ind3   = 5;
           ind4   = 6;
           ind5   = 7;
           ind6   = 8;
           n = alpha + beta * FH;
        else
            warning('Negative number of feedheaters not allowed')
        end
    case 'on'
        if FH == 0
            alpha = 8;
            ind3HP  = 5;
            ind4HP  = 6;
            ind3LP  = 7;
            ind4LP  = 8;
        elseif FH > 0
            alpha = 12;
            beta  = 4;
            ind3HP  = 5;
            ind4HP  = 6;
            ind3LP  = 7;
            ind4LP  = 8;
            ind5    = 9;
            ind6    = 10;
            n = alpha + beta * FH + 4; 
        else
            warning('Negative number of feedheaters not allowed')
        end
    otherwise
        warning('Unexpected reheating entry.')
end


% Determining the Bache Position
if FH > 4
    index_bleed_bache = alpha + 1;
    position_bache = 1;
    while state{index_bleed_bache}.p < 4.6
        position_bache = position_bache + 1;
        index_bleed_bache = index_bleed_bache + beta;
    end
end

%% Saturated States Curve
S = linspace(0.00139,9.1531,2*step2);
Tsat = arrayfun(@(s) XSteam('Tsat_s',s), S);
Psat = arrayfun(@(s) XSteam('psat_s',s), S); 
Hsat = arrayfun(@(p,s) XSteam('h_ps',p,s), Psat, S); % entalpie saturee pour chaque entropie

%% Main Commom Points
% Curve 1 to 2
s12 = linspace(state{ind1}.s,state{ind2}.s,step2);
p12 = linspace(state{ind1}.p,state{ind2}.p,step2);
h12 = arrayfun(@(p,s) XSteam('h_ps',p,s), p12, s12);
t12 = arrayfun(@(p,s) XSteam('T_ps',p,s), p12, s12);
% Curve 2 to 2i
s22i = linspace(state{ind2}.s,state{ind2i}.s,step2);
p22i = linspace(state{ind2}.p,state{ind2i}.p,step2);
h22i = arrayfun(@(p,s) XSteam('h_ps',p,s), p22i, s22i);
t22i = arrayfun(@(p,s) XSteam('T_ps',p,s), p22i, s22i);
% Curve 2i to 2ii
s2i2ii = linspace(state{ind2i}.s,state{ind2ii}.s,step2);
p2i2ii = linspace(state{ind2i}.p,state{ind2ii}.p,step2);
h2i2ii = arrayfun(@(p,s) XSteam('h_ps',p,s), p2i2ii, s2i2ii);
t2i2ii = arrayfun(@(p,s) XSteam('T_ps',p,s), p2i2ii, s2i2ii);


s12ii = [s12, s22i, s2i2ii];
t12ii = [t12, t22i, t2i2ii];
h12ii = [h12, h22i, h2i2ii];

%% RH OFF, FH = 0
if alpha == 6
    % Curve 2'' to 3
    s2ii3 = linspace(state{ind2ii}.s,state{ind3}.s,step2);
    p2ii3 = linspace(state{ind2ii}.p,state{ind3}.p,step2);
    h2ii3 = arrayfun(@(p,s) XSteam('h_ps',p,s), p2ii3, s2ii3);
    t2ii3 = arrayfun(@(p,s) XSteam('T_ps',p,s), p2ii3, s2ii3);
    % Curve 3 to 4
    p34   = linspace(state{ind3}.p,state{ind4}.p,step2);
    s34_s = linspace(state{ind3}.s,state{ind3}.s,step2);
    h34_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34,s34_s);
    h34   = zeros(1,step2);
    for i = 1:step2
        h34(i) = state{ind3}.h + eta_SiTLP*(h34_s(i) - state{ind3}.h);
    end
    t34   = arrayfun(@(p,h) XSteam('T_ph',p,h),p34,h34);
    s34   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34,h34);

    % Curve Inside Condensor 
    % 4 to 4L (vapor to sat liquid)
    s4L  = XSteam('sL_p',state{ind4}.p);
    s44L = linspace(state{ind4}.s,s4L,step2);
    p44L = linspace(state{ind4}.p,state{ind4}.p,step2);
    h44L = arrayfun(@(p,s) XSteam('h_ps',p,s), p44L, s44L);
    t44L = arrayfun(@(p,s) XSteam('T_ps',p,s), p44L, s44L);
    % 4L to 1 sat liq to undercooled liquid at exit of condensor
    s4L1 = linspace(s4L,state{ind1}.s,step2);
    p4L1 = linspace(state{ind4}.p,state{ind1}.p,step2);
    h4L1 = arrayfun(@(p,s) XSteam('h_ps',p,s), p4L1, s4L1);
    t4L1 = arrayfun(@(p,s) XSteam('T_ps',p,s), p4L1, s4L1);
    
    s41 = [s44L, s4L1];
    t41 = [t44L, t4L1];
    h41 = [h44L, h4L1];
    
Sif = [s2ii3, s34, s41];
Tif = [t2ii3, t34, t41];
Hif = [h2ii3, h34, h41];

%% RH ON, FH = 0
elseif alpha == 8 
    % Main Points
    % Curve 2'' to 3HP
    s2ii3 = linspace(state{ind2ii}.s,state{ind3HP}.s,step2);
    p2ii3 = linspace(state{ind2ii}.p,state{ind3HP}.p,step2);
    h2ii3 = arrayfun(@(p,s) XSteam('h_ps',p,s), p2ii3, s2ii3);
    t2ii3 = arrayfun(@(p,s) XSteam('T_ps',p,s), p2ii3, s2ii3);
    
    % Curve 3HP to 4HP
    p34HP   = linspace(state{ind3HP}.p,state{ind4HP}.p,step2);
    s34HP_s = linspace(state{ind3HP}.s,state{ind3HP}.s,step2);
    h34HP_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34HP,s34HP_s);
    h34HP   = zeros(1,step2);
    for i = 1:step2
        h34HP(i) = state{ind3HP}.h + eta_SiTHP*(h34HP_s(i) - state{ind3HP}.h);
    end
    t34HP   = arrayfun(@(p,h) XSteam('T_ph',p,h),p34HP,h34HP);
    s34HP   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34HP,h34HP);
    % Curve 4HP to 3LP
    s4HP3LP = linspace(state{ind4HP}.s,state{ind3LP}.s,step2);
    p4HP3LP = linspace(state{ind4HP}.p,state{ind3LP}.p,step2);
    h4HP3LP = arrayfun(@(p,s) XSteam('h_ps',p,s), p4HP3LP, s4HP3LP);
    t4HP3LP = arrayfun(@(p,s) XSteam('T_ps',p,s), p4HP3LP, s4HP3LP);  
    % Curve 3BP to 4BP
    p34LP   = linspace(state{ind3LP}.p,state{ind4LP}.p,step2);
    s34LP_s = linspace(state{ind3LP}.s,state{ind3LP}.s,step2);
    h34LP_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34LP,s34LP_s);
    h34LP   = zeros(1,step2);
    for i = 1:step2
        h34LP(i) = state{ind3LP}.h + eta_SiTLP*(h34LP_s(i) - state{ind3LP}.h);
    end
    t34LP   = arrayfun(@(p,h) XSteam('T_ph',p,h),p34LP,h34LP);
    s34LP   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34LP,h34LP);
    
    s34 = [s34HP, s4HP3LP, s34LP];
    t34 = [t34HP, t4HP3LP, t34LP];
    h34 = [h34HP, h4HP3LP, h34LP];

    % Curve Inside Condensor 
    % 4LP to 4L (vapor to sat liquid)
    s4LPL  = XSteam('sL_p',state{ind4LP}.p);
    s4LP4L = linspace(state{ind4LP}.s,s4LPL,step2);
    p4LP4L = linspace(state{ind4LP}.p,state{ind4LP}.p,step2);
    h4LP4L = arrayfun(@(p,s) XSteam('h_ps',p,s), p4LP4L, s4LP4L);
    t4LP4L = arrayfun(@(p,s) XSteam('T_ps',p,s), p4LP4L, s4LP4L);
    % 4L to 1 sat liq to undercooled liquid at exit of condensor
    s4L1 = linspace(s4LPL,state{ind1}.s,step2);
    p4L1 = linspace(state{ind4LP}.p,state{ind1}.p,step2);
    h4L1 = arrayfun(@(p,s) XSteam('h_ps',p,s), p4L1, s4L1);
    t4L1 = arrayfun(@(p,s) XSteam('T_ps',p,s), p4L1, s4L1);
    
    s4LP1 = [s4LP4L, s4L1];
    t4LP1 = [t4LP4L, t4L1];
    h4LP1 = [h4LP4L, h4L1];
    
    
    Sif = [s2ii3, s34, s4LP1];
    Tif = [t2ii3, t34, t4LP1];
    Hif = [h2ii3, h34, h4LP1];
    
%% RH OFF, FH # 0
elseif alpha == 10
    % Main Points
    % Curve 2'' to 3
    s2ii3 = linspace(state{ind2ii}.s,state{ind3}.s,step2);
    p2ii3 = linspace(state{ind2ii}.p,state{ind3}.p,step2);
    h2ii3 = arrayfun(@(p,s) XSteam('h_ps',p,s), p2ii3, s2ii3);
    t2ii3 = arrayfun(@(p,s) XSteam('T_ps',p,s), p2ii3, s2ii3);
    % Curve 3 to 4
    p34   = linspace(state{ind3}.p,state{ind4}.p,step2);
    s34_s = linspace(state{ind3}.s,state{ind3}.s,step2);
    h34_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34,s34_s);
    h34   = zeros(1,step2);
    for i = 1:step2
        h34(i) = state{ind3}.h + eta_SiTLP*(h34_s(i) - state{ind3}.h);
    end
    t34   = arrayfun(@(p,h) XSteam('T_ph',p,h),p34,h34);
    s34   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34,h34);
    
    % Curve Inside Condensor 
    % 4 to 4L (vapor to sat liquid)
    s4L = XSteam('sL_p',state{ind4}.p);
    s44L = linspace(state{ind4}.s,s4L,step2);
    p44L = linspace(state{ind4}.p,state{ind4}.p,step2);
    h44L = arrayfun(@(p,s) XSteam('h_ps',p,s), p44L, s44L);
    t44L = arrayfun(@(p,s) XSteam('T_ps',p,s), p44L, s44L);
    % 4L to 5 sat liq to undercooled liquid at exit of condensor
    s4L5 = linspace(s4L,state{ind5}.s,step2);
    p4L5 = linspace(state{ind4}.p,state{ind5}.p,step2);
    h4L5 = arrayfun(@(p,s) XSteam('h_ps',p,s), p4L5, s4L5);
    t4L5 = arrayfun(@(p,s) XSteam('T_ps',p,s), p4L5, s4L5);
    
    s45 = [s44L, s4L5];
    t45 = [t44L, t4L5];
    h45 = [h44L, h4L5];
    
    % Condensor Pump 5 to 6
    s56 = linspace(state{ind5}.s,state{ind6}.s,step2);
    p56 = linspace(state{ind5}.p,state{ind6}.p,step2);
    h56 = arrayfun(@(p,s) XSteam('h_ps',p,s), p56, s56);
    t56 = arrayfun(@(p,s) XSteam('T_ps',p,s), p56, s56);
    % Curves from 6 to 1
    if FH > 4
        index_bleed_bache = alpha + 1;
        position_bache = 1;
        while state{index_bleed_bache}.p < 4.6
            position_bache = position_bache + 1;
            index_bleed_bache = index_bleed_bache + beta;
        end
        % from 6 to bache
        s6bache = linspace(state{ind6}.s,state{index_bleed_bache + 2}.s,step2);
        p6bache = linspace(state{ind6}.p,state{index_bleed_bache + 2}.p,step2);
        h6bache = arrayfun(@(p,s) XSteam('h_ps',p,s), p6bache, s6bache);
        t6bache = arrayfun(@(p,s) XSteam('T_ps',p,s), p6bache, s6bache);
        
        % bache pump
        ind6BP = alpha + beta*(position_bache + 1) - 2;
        sBP = linspace(state{index_bleed_bache + 2}.s,state{ind6BP}.s,step2);
        pBP = linspace(state{index_bleed_bache + 2}.p,state{ind6BP}.p,step2);
        hBP = arrayfun(@(p,s) XSteam('h_ps',p,s), pBP, sBP);
        tBP = arrayfun(@(p,s) XSteam('T_ps',p,s), pBP, sBP);
        
        % from 6BP to 1
        s6BP1 = linspace(state{ind6BP}.s,state{ind1}.s,step2);
        p6BP1 = linspace(state{ind6BP}.p,state{ind1}.p,step2);
        h6BP1 = arrayfun(@(p,s) XSteam('h_ps',p,s), p6BP1, s6BP1);
        t6BP1 = arrayfun(@(p,s) XSteam('T_ps',p,s), p6BP1, s6BP1);
        
        s61 = [s6bache, sBP, s6BP1];
        t61 = [t6bache, tBP, t6BP1];
        h61 = [h6bache, hBP, h6BP1];
    else
        s61 = linspace(state{ind6}.s,state{ind1}.s,step2);
        p61 = linspace(state{ind6}.p,state{ind1}.p,step2);
        t61 = arrayfun(@(p,s) XSteam('T_ps',p,s), p61, s61); 
        h61 = arrayfun(@(p,s) XSteam('h_ps',p,s), p61, s61);
    end

    Sif = [s2ii3, s34, s45, s56, s61];
    Tif = [t2ii3, t34, t45, t56, t61];
    Hif = [h2ii3, h34, h45, h56, h61];
    
    % Additional Points
    % Curves 4i to 7i 
    s47i = zeros(FH,step2);
    p47i = zeros(FH,step2);
    h47i = zeros(FH,step2);
    t47i = zeros(FH,step2);
    for i = 1:FH
        ind4i = alpha + beta*(i-1) + 1;
        ind7i = alpha + beta*i - 1;
        if (FH > 4) && (i == position_bache)
            p_bleedbache = state{ind7i}.p; % Pressure behind the isenthalpic valve 
            h_bleedbache = state{ind4i}.h; 
            s_bleedbache = XSteam('s_ph',p_bleedbache,h_bleedbache);
            s47i(i,:) = linspace(s_bleedbache,state{ind7i}.s,step2);
            p47i(i,:) = linspace(p_bleedbache,state{ind7i}.p,step2);
            h47i(i,:) = arrayfun(@(p,s) XSteam('h_ps',p,s), p47i(i,:), s47i(i,:));
            t47i(i,:) = arrayfun(@(p,s) XSteam('T_ps',p,s), p47i(i,:), s47i(i,:)); 
            
            s4valve = linspace(state{ind4i}.s,s_bleedbache,step2);
            p4valve = linspace(state{ind4i}.p,p_bleedbache,step2);
            h4valve = arrayfun(@(p,s) XSteam('h_ps',p,s), p4valve, s4valve);
            t4valve = arrayfun(@(p,s) XSteam('T_ps',p,s), p4valve, s4valve);
        else
        s47i(i,:) = linspace(state{ind4i}.s,state{ind7i}.s,step2);
        p47i(i,:) = linspace(state{ind4i}.p,state{ind7i}.p,step2);
        h47i(i,:) = arrayfun(@(p,s) XSteam('h_ps',p,s), p47i(i,:), s47i(i,:));
        t47i(i,:) = arrayfun(@(p,s) XSteam('T_ps',p,s), p47i(i,:), s47i(i,:)); 
        end
    end
    % Additional curves were intentionally left out because they would not
    % bring much interesting information and would render the calculations
    % heavier.
    
    
%% RH ON, FH # 0
elseif alpha == 12
    % Main Points
    % Curve 2'' to 3HP
    s2ii3 = linspace(state{ind2ii}.s,state{ind3HP}.s,step2);
    p2ii3 = linspace(state{ind2ii}.p,state{ind3HP}.p,step2);
    h2ii3 = arrayfun(@(p,s) XSteam('h_ps',p,s), p2ii3, s2ii3);
    t2ii3 = arrayfun(@(p,s) XSteam('T_ps',p,s), p2ii3, s2ii3);
    
    % Curve 3HP to 4HP
    p34HP   = linspace(state{ind3HP}.p,state{ind4HP}.p,step2);
    s34HP_s = linspace(state{ind3HP}.s,state{ind3HP}.s,step2);
    h34HP_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34HP,s34HP_s);
    h34HP   = zeros(1,step2);
    for i = 1:step2
        h34HP(i) = state{ind3HP}.h + eta_SiTHP*(h34HP_s(i) - state{ind3HP}.h);
    end
    t34HP   = arrayfun(@(p,h) XSteam('T_ph',p,h),p34HP,h34HP);
    s34HP   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34HP,h34HP);
    % Curve 4HP to 3LP
    s4HP3LP = linspace(state{ind4HP}.s,state{ind3LP}.s,step2);
    p4HP3LP = linspace(state{ind4HP}.p,state{ind3LP}.p,step2);
    h4HP3LP = arrayfun(@(p,s) XSteam('h_ps',p,s), p4HP3LP, s4HP3LP);
    t4HP3LP = arrayfun(@(p,s) XSteam('T_ps',p,s), p4HP3LP, s4HP3LP);  
    % Curve 3BP to 4BP
    p34LP   = linspace(state{ind3LP}.p,state{ind4LP}.p,step2);
    s34LP_s = linspace(state{ind3LP}.s,state{ind3LP}.s,step2);
    h34LP_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34LP,s34LP_s);
    h34LP   = zeros(1,step2);
    for i = 1:step2
        h34LP(i) = state{ind3LP}.h + eta_SiTLP*(h34LP_s(i) - state{ind3LP}.h);
    end
    t34LP   = arrayfun(@(p,h) XSteam('T_ph',p,h),p34LP,h34LP);
    s34LP   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34LP,h34LP);
        
    s34 = [s34HP, s4HP3LP, s34LP];
    t34 = [t34HP, t4HP3LP, t34LP];
    h34 = [h34HP, h4HP3LP, h34LP];
    
    % Curve Inside Condensor 
    % 4LP to 4L (vapor to sat liquid)
    s4LPL  = XSteam('sL_p',state{ind4LP}.p);
    s4LP4L = linspace(state{ind4LP}.s,s4LPL,step2);
    p4LP4L = linspace(state{ind4LP}.p,state{ind4LP}.p,step2);
    h4LP4L = arrayfun(@(p,s) XSteam('h_ps',p,s), p4LP4L, s4LP4L);
    t4LP4L = arrayfun(@(p,s) XSteam('T_ps',p,s), p4LP4L, s4LP4L);
    % 4L to 1 sat liq to undercooled liquid at exit of condensor
    s4L5 = linspace(s4LPL,state{ind5}.s,step2);
    p4L5 = linspace(state{ind4LP}.p,state{ind5}.p,step2);
    h4L5 = arrayfun(@(p,s) XSteam('h_ps',p,s), p4L5, s4L5);
    t4L5 = arrayfun(@(p,s) XSteam('T_ps',p,s), p4L5, s4L5);
    
    s45 = [s4LP4L, s4L5];
    t45 = [t4LP4L, t4L5];
    h45 = [h4LP4L, h4L5];
    
    % Condensor Pump 1
    s56 = linspace(state{ind5}.s,state{ind6}.s,step2);
    p56 = linspace(state{ind5}.p,state{ind6}.p,step2);
    h56 = arrayfun(@(p,s) XSteam('h_ps',p,s), p56, s56);
    t56 = arrayfun(@(p,s) XSteam('T_ps',p,s), p56, s56);
    % Curves from 6 to 1
    if FH > 4
        % from 6 to bache
        s6bache = linspace(state{ind6}.s,state{index_bleed_bache + 2}.s,step2);
        p6bache = linspace(state{ind6}.p,state{index_bleed_bache + 2}.p,step2);
        h6bache = arrayfun(@(p,s) XSteam('h_ps',p,s), p6bache, s6bache);
        t6bache = arrayfun(@(p,s) XSteam('T_ps',p,s), p6bache, s6bache);
        
        % bache pump
        ind6BP = alpha + beta*(position_bache + 1) - 2;
        sBP = linspace(state{index_bleed_bache + 2}.s,state{ind6BP}.s,step2);
        pBP = linspace(state{index_bleed_bache + 2}.p,state{ind6BP}.p,step2);
        hBP = arrayfun(@(p,s) XSteam('h_ps',p,s), pBP, sBP);
        tBP = arrayfun(@(p,s) XSteam('T_ps',p,s), pBP, sBP);
        
        % from 6BP to 1
        s6BP1 = linspace(state{ind6BP}.s,state{ind1}.s,step2);
        p6BP1 = linspace(state{ind6BP}.p,state{ind1}.p,step2);
        h6BP1 = arrayfun(@(p,s) XSteam('h_ps',p,s), p6BP1, s6BP1);
        t6BP1 = arrayfun(@(p,s) XSteam('T_ps',p,s), p6BP1, s6BP1);
        
        s61 = [s6bache, sBP, s6BP1];
        t61 = [t6bache, tBP, t6BP1];
        h61 = [h6bache, hBP, h6BP1];
    else
        s61 = linspace(state{ind6}.s,state{ind1}.s,step2);
        p61 = linspace(state{ind6}.p,state{ind1}.p,step2);
        t61 = arrayfun(@(p,s) XSteam('T_ps',p,s), p61, s61); 
        h61 = arrayfun(@(p,s) XSteam('h_ps',p,s), p61, s61);
    end
    
    Sif = [s2ii3, s34, s45, s56, s61];
    Tif = [t2ii3, t34, t45, t56, t61];
    Hif = [h2ii3, h34, h45, h56, h61];
    
    % Additional Points
    % Curve 4HP to 7HP
    s4HP7 = linspace(state{ind4HP}.s,state{n - beta + 3}.s,step2);
    p4HP7 = linspace(state{ind4HP}.p,state{n - beta + 3}.p,step2);
    h4HP7 = arrayfun(@(p,s) XSteam('h_ps',p,s), p4HP7, s4HP7);
    t4HP7 = arrayfun(@(p,s) XSteam('T_ps',p,s), p4HP7, s4HP7);

    % Curves 4i to 7i 
    s47i = zeros(FH + 1,step2);
    p47i = zeros(FH + 1,step2);
    h47i = zeros(FH + 1,step2);
    t47i = zeros(FH + 1,step2);
    for i = 1:FH
        ind4i = alpha + beta*(i-1) + 1;
        ind7i = alpha + beta*i - 1;
        if (FH > 4) && (i == position_bache)
            p_bleedbache = state{ind7i}.p; % Pressure behind the isenthalpic valve 
            h_bleedbache = state{ind4i}.h; 
            s_bleedbache = XSteam('s_ph',p_bleedbache,h_bleedbache);
            s47i(i,:) = linspace(s_bleedbache,state{ind7i}.s,step2);
            p47i(i,:) = linspace(p_bleedbache,state{ind7i}.p,step2);
            h47i(i,:) = arrayfun(@(p,s) XSteam('h_ps',p,s), p47i(i,:), s47i(i,:));
            t47i(i,:) = arrayfun(@(p,s) XSteam('T_ps',p,s), p47i(i,:), s47i(i,:)); 
            
            s4valve = linspace(state{ind4i}.s,s_bleedbache,step2);
            p4valve = linspace(state{ind4i}.p,p_bleedbache,step2);
            h4valve = arrayfun(@(p,s) XSteam('h_ps',p,s), p4valve, s4valve);
            t4valve = arrayfun(@(p,s) XSteam('T_ps',p,s), p4valve, s4valve);
        else
        s47i(i,:) = linspace(state{ind4i}.s,state{ind7i}.s,step2);
        p47i(i,:) = linspace(state{ind4i}.p,state{ind7i}.p,step2);
        h47i(i,:) = arrayfun(@(p,s) XSteam('h_ps',p,s), p47i(i,:), s47i(i,:));
        t47i(i,:) = arrayfun(@(p,s) XSteam('T_ps',p,s), p47i(i,:), s47i(i,:)); 
        end
    end
    s47i(FH+1,:) = s4HP7;
    t47i(FH+1,:) = t4HP7;
    h47i(FH+1,:) = h4HP7;
    % Additional curves were intentionally left out because they would not
    % bring much interesting information and would render the calculations
    % heavier.
    
end


%% Individual Cycle Points
A = size(state);
if FH > 0
    T_Cycle = ones(1,A(1) - 1);
    H_Cycle = ones(1,A(1) - 1);
    S_Cycle = ones(1,A(1) - 1);
    for i=1:A(1) - 1
        if i < alpha + 2
            T_Cycle(i) = state{i}.t;
            H_Cycle(i) = state{i}.h;
            S_Cycle(i) = state{i}.s;
        elseif i >= alpha + 2
            T_Cycle(i) = state{i + 1}.t;
            H_Cycle(i) = state{i + 1}.h;
            S_Cycle(i) = state{i + 1}.s;
        end
    end
else
    T_Cycle = ones(1,A(1));
    H_Cycle = ones(1,A(1));
    S_Cycle = ones(1,A(1));
    for i = 1:A(1)
        T_Cycle(i) = state{i}.t;
        H_Cycle(i) = state{i}.h;
        S_Cycle(i) = state{i}.s;
    end
end



Smain = [s12ii, Sif];
Tmain = [t12ii, Tif];
Hmain = [h12ii, Hif];


figure
subplot(2,1,1);
plot(S_Cycle,T_Cycle,'*r',S,Tsat,'-b',Smain,Tmain,'-r'); hold on;
if FH > 0
    size_s47i = size(s47i);
    for i = 1:size_s47i(1)
        plot(s47i(i,:),t47i(i,:),'-.g'); hold on;
    end
    if FH > 4
        plot(s4valve,t4valve,'-.r'); hold on;
    end
end
xlabel('s [kJ/kgK]'); hold on;
ylabel('t [°C]'); hold on;
grid on; 

subplot(2,1,2);
plot(S_Cycle,H_Cycle,'*r',S,Hsat,'-b',Smain,Hmain,'-r'); hold on;
if FH > 0
    size_s47i = size(s47i);
    for i = 1:size_s47i(1)
        plot(s47i(i,:),h47i(i,:),'-.g'); hold on;
    end
    if FH > 4
        plot(s4valve,h4valve,'-.r'); hold on;
    end
end
xlabel('s [kJ/kgK]'); hold on;
ylabel('h [kJ/kg]'); hold on;
grid on; 
end

