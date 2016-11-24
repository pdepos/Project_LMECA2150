function [ Output ] = RankineHirnOLD( t3, p3 )
%Rankine-Hirn: Simple Rankine-Hirn cycle.
%   Input Arguments: 
%   - Temperature at point 3 [°C]
%   - Pressure at point 3 [bar]
%   =======================================================================
t4 = 32;
t1 = t4;

%Reference State
t0 = 273.15 + 15; %[K]
p0 = 0.01704;
h0 = 63.0;
s0 = 0.224;

p = [1;1;p3;1];
t = [t1;1;t3;t4];
x = [0;NaN;NaN;1];
h = ones(4,1);
s = ones(4,1);
e = ones(4,1);
v_LH2O   = 0.001005; %[m³/kg] volume massique de l'eau
eta_SiT = 0.88; % Isentropic Efficiency turbine
eta_iP  = 0.85; % Efficiency feedwater pump
kpdgen  = 0.90; % Pressure drop coefficient at steam generator. Determined via numerical examples in book


% Point 3
h(3) = XSteam('h_pT',p(3),t(3));
s(3) = XSteam ('s_pT',p(3),t(3));
e(3) = h(3) - h0 - t0*(s(3) - s0);

% Point 4
p(4) = XSteam('psat_t',t(4));
s4s  = s(3);
x4s  = (s4s - XSteam('sL_T',t(4))) / (XSteam('sV_T',t(4)) - XSteam('sL_T',t(4)));
h4s  = x4s*XSteam('hV_T',t(4)) + (1-x4s)*XSteam('hL_T',t(4));
h(4) = eta_SiT*(h4s - h(3)) + h(3);
x(4) = (h(4) - XSteam('hL_T',t(4))) / (XSteam('hV_T',t(4)) - XSteam('hL_T',t(4)));
s(4) = x(4)*XSteam('sV_T',t(4)) + (1-x(4))*XSteam('sL_T',t(4));
e(4) = h(4) - h0 - t0*(s(4) - s0);

% Point 1
p(1) = p(4);
h(1) = XSteam('hL_T',t(1));
s(1) = XSteam('sL_t',t(1));
e(1) = h(1) - h0 - t0*(s(1) - s0);

% Point 2
p(2) = p(3)/kpdgen;
h(2) = v_LH2O * (p(2) - p(1))*1E2 / eta_iP + h(1); % *e2 pour obtenir kJ/kg (si e5 on obtient des joules..)
%muT2 = (XSteam('my_pT',p(1),t(1)) + XSteam('my_pT',p(2),t(1))) / 2 %muT première approximation
muT2 = 0.0905;
cp2  = (XSteam('Cp_pT',p(2),t(1)) + XSteam('CpL_p',p(1))) / 2; %Output 
t(2) =  t(1) + (v_LH2O*1e2 - muT2)*(p(2) - p(1))/ (cp2);
% Should add iterations to refine point 2
s(2) = XSteam ('s_pT',p(2),t(2));
e(2) = h(2) - h0 - t0*(s(2) - s0);


Output = [t p x h s e];
end

