function [sL,sV,hL,hV,s115,H115,s318,H318,s818,H818,s1920,h1920,s201,h201] = combined3hs(stateV)


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
[s12,h12] = transitionHS(stateV,1,2,'undefined',500);

%%%% 2 --> 3 (isobaric) %%%%
[s23,h23] = transitionHS(stateV,2,3,'isobaric',500);

%%%% 3 --> 4 (saturation) %%%%
[s34,h34] = transitionHS(stateV,3,4,'saturation',500);

%%%% 4 --> 5 (isobaric) %%%%
[s45,h45] = transitionHS(stateV,4,5,'isobaric',500);

%%%% 5 --> 6 (isobaric) %%%%
[s56,h56] = transitionHS(stateV,5,6,'isobaric',500);

%%%% 3 --> 7 (pump) %%%%
[s37,h37] = transitionHS(stateV,3,7,'undefined',500);

%%%% 7 --> 8 (isobaric) %%%%
[s78,h78] = transitionHS(stateV,7,8,'isobaric',500);

%%%% 8 --> 12 (saturation) %%%%
[s812,h812] = transitionHS(stateV,8,12,'saturation',500);

%%%% 12 --> 13 (isobaric) %%%%
[s1213,h1213] = transitionHS(stateV,12,13,'isobaric',500);

%%%% 8 --> 9 (pump) %%%%
[s89,h89] = transitionHS(stateV,8,9,'undefined',500);

%%%% 9 --> 10 (isobaric) %%%%
[s910,h910] = transitionHS(stateV,9,10,'isobaric',500);

%%%% 10 --> 11 (saturation) %%%%
[s1011,h1011] = transitionHS(stateV,10,11,'saturation',500);

%%%% 11 --> 14 (isobaric) %%%%
[s1114,h1114] = transitionHS(stateV,11,14,'isobaric',500);

%%%% 14 --> 15 (turbine) %%%%
[s1415,h1415] = transitionHS(stateV,14,15,'undefined',500);

%%%% 16 --> 17 (isobaric) %%%%
[s1617,h1617] = transitionHS(stateV,16,17,'isobaric',500);

%%%% 17 --> 18 (turbine) %%%%
[s1718,h1718] = transitionHS(stateV,17,18,'undefined',500);

%%%% 19 --> 20 (turbine) %%%%
[s1920,h1920] = transitionHS(stateV,19,20,'undefined',500);

%%%% 6 --> 18 (isobaric) %%%%
[s618,h618] = transitionHS(stateV,6,18,'isobaric',500);

%%%% 20 --> 1 (saturation) %%%% 
[s201,h201] = transitionHS(stateV,1,20,'saturation',500);

%%%%%%%%%%%%%%%%%%
%%%% Assembly %%%%
%%%%%%%%%%%%%%%%%%

s115 = [s12 s23 s37 s78 s89 s910 s1011 s1114 s1415];
s318 = [s34 s45 s56 s618];
s818 = [s812 s1213 s1617 s1718];
%s191 = [s201 s1920];

H115 = [h12 h23 h37 h78 h89 h910 h1011 h1114 h1415];
H318 = [h34 h45 h56 h618];
H818 = [h812 h1213 h1617 h1718];
%H191 = [h201 h1920];

end