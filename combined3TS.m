function [sL,sV,TL,TV,s115,T115,s318,T318,s818,T818,s191,T191] = combined3TS(stateV)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot the saturation curve %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ts = linspace(0,1000,100000);
lt = length(Ts);

sL = zeros(lt,1);
sV = zeros(lt,1);

for i = 1:lt
   sL(i) = XSteam('sL_T',Ts(i));
   sV(i) = XSteam('sV_T',Ts(i));
end

indL = find(isfinite(sL));
indV = find(isfinite(sV));

sL = sL(indL);
sV = sV(indV);

TL = Ts(indL);
TV = Ts(indV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot the transitions %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 1 --> 2 %%%%
[s12,T12] = transitionTS(stateV,1,2,'undefined',500);

%%%% 2 --> 3 (isobaric) %%%%
[s23,T23] = transitionTS(stateV,2,3,'isobaric',500);

%%%% 3 --> 4 (saturation) %%%%
[s34,T34] = transitionTS(stateV,3,4,'saturation',500);

%%%% 4 --> 5 (isobaric) %%%%
[s45,T45] = transitionTS(stateV,4,5,'isobaric',500);

%%%% 5 --> 6 (isobaric) %%%%
[s56,T56] = transitionTS(stateV,5,6,'isobaric',500);

%%%% 3 --> 7 (pump) %%%%
[s37,T37] = transitionTS(stateV,3,7,'undefined',500);

%%%% 7 --> 8 (isobaric) %%%%
[s78,T78] = transitionTS(stateV,7,8,'isobaric',500);

%%%% 8 --> 12 (saturation) %%%%
[s812,T812] = transitionTS(stateV,8,12,'saturation',500);

%%%% 12 --> 13 (isobaric) %%%%
[s1213,T1213] = transitionTS(stateV,12,13,'isobaric',500);

%%%% 8 --> 9 (pump) %%%%
[s89,T89] = transitionTS(stateV,8,9,'undefined',500);

%%%% 9 --> 10 (isobaric) %%%%
[s910,T910] = transitionTS(stateV,9,10,'isobaric',500);

%%%% 10 --> 11 (saturation) %%%%
[s1011,T1011] = transitionTS(stateV,10,11,'saturation',500);

%%%% 11 --> 14 (isobaric) %%%%
[s1114,T1114] = transitionTS(stateV,11,14,'isobaric',500);

%%%% 14 --> 15 (turbine) %%%%
[s1415,T1415] = transitionTS(stateV,14,15,'undefined',500);

%%%% 16 --> 17 (isobaric) %%%%
[s1617,T1617] = transitionTS(stateV,16,17,'isobaric',500);

%%%% 17 --> 18 (turbine) %%%%
[s1718,T1718] = transitionTS(stateV,17,18,'undefined',500);

%%%% 19 --> 20 (turbine) %%%%
[s1920,T1920] = transitionTS(stateV,19,20,'undefined',500);

%%%% 6 --> 18 (isobaric) %%%%
[s618,T618] = transitionTS(stateV,6,18,'isobaric',500);

%%%% 20 --> 1 (saturation) %%%% 
[s201,T201] = transitionTS(stateV,20,1,'saturation',500);

%%%%%%%%%%%%%%%%%%
%%%% Assembly %%%%
%%%%%%%%%%%%%%%%%%

s115 = [s12 s23 s37 s78 s89 s910 s1011 s1114 s1415];
s318 = [s34 s45 s56 s618];
s818 = [s812 s1213 s1617 s1718];
s191 = [s1920 s201];

T115 = [T12 T23 T37 T78 T89 T910 T1011 T1114 T1415];
T318 = [T34 T45 T56 T618];
T818 = [T812 T1213 T1617 T1718];
T191 = [T1920 T201];

end