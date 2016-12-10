function [sL,sV,hL,hV,s110,h110,s39,h39] = combined2hs(stateV)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot the saturation curve %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = linspace(0,5000,25000);
lT = length(T);

sL = zeros(lT,1);
sV = zeros(lT,1);
hL = zeros(lT,1);
hV = zeros(lT,1);

for i = 1:lT
   sL(i) = XSteam('sL_T',T(i));
   sV(i) = XSteam('sV_T',T(i));
   hL(i) = XSteam('hL_T',T(i));
   hV(i) = XSteam('hV_T',T(i));
end

indsL = find(isfinite(sL));
indsV = find(isfinite(sV));
indhL = find(isfinite(hL));
indhV = find(isfinite(hV));

sL = sL(indsL);
sV = sV(indsV);
hL = hL(indhL);
hV = hV(indhV);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot transitions %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 1 --> 2 %%%%
[s12,h12] = transitionHS(stateV,1,2,'undefined',200);

%%%% 2 --> 3 (isobaric) %%%%
[s23,h23] = transitionHS(stateV,2,3,'isobaric',300);

%%%% 3 --> 4 (isobaric and iso-T) %%%%
[s34,h34] = transitionHS(stateV,3,4,'saturation',500);


%%%% 4 --> 9 (isobaric) %%%%
[s49,h49] = transitionHS(stateV,4,9,'isobaric',500);

%%%% 9 --> 10 %%%%
[s910,h910] = transitionHS(stateV,9,10,'undefined',1000);

%%%% 3 --> 5 %%%%
[s35,h35] = transitionHS(stateV,3,5,'undefined',200);

%%%% 5 --> 6 (isobaric) %%%%
[s56,h56] = transitionHS(stateV,5,6,'isobaric',500);

%%%% 6 --> 7 (saturation) %%%%
[s67,h67] = transitionHS(stateV,6,7,'saturation',500);

%%%% 7 --> 8 (isobaric) %%%%
[s78,h78] = transitionHS(stateV,7,8,'isobaric',1000);

%%%% 8 --> 9 %%%%
[s89,h89] = transitionHS(stateV,8,9,'undefined',1000);

%%%% 10 --> 1 %%%%
[s101,h101] = transitionHS(stateV,10,1,'saturation',500);

%%%% Assembly %%%%
s110 = [s12 s23 s34 s49 s910 s101];
s39 = [s35 s56 s67 s78 s89];

h110 = [h12 h23 h34 h49 h910 h101];
h39 = [h35 h56 h67 h78 h89];

end