function [sL,sV,TL,TV,s110,T110,s39,T39] = combined2TS(stateV)

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
[s12,T12] = transitionTS(stateV,1,2,'undefined',200);

%%%% 2 --> 3 (isobaric) %%%%
[s23,T23] = transitionTS(stateV,2,3,'isobaric',300);

%%%% 3 --> 4 (isobaric and iso-T) %%%%
[s34,T34] = transitionTS(stateV,3,4,'saturation',500);

%%%% 4 --> 9 (isobaric) %%%%
[s49,T49] = transitionTS(stateV,4,9,'isobaric',500);

%%%% 9 --> 10 %%%%
[s910,T910] = transitionTS(stateV,9,10,'undefined',1000);

%%%% 3 --> 5 %%%%
[s35,T35] = transitionTS(stateV,3,5,'undefined',200);

%%%% 5 --> 6 (isobaric) %%%%
[s56,T56] = transitionTS(stateV,5,6,'isobaric',500);

%%%% 6 --> 7 (saturation) %%%%
[s67,T67] = transitionTS(stateV,6,7,'saturation',500);

%%%% 7 --> 8 (isobaric) %%%%
[s78,T78] = transitionTS(stateV,7,8,'isobaric',1000);

%%%% 8 --> 9 %%%%
[s89,T89] = transitionTS(stateV,8,9,'undefined',1000);

%%%% 10 --> 1 %%%%
[s101,T101] = transitionTS(stateV,10,1,'saturation',500);

%%%% Assembly %%%%
s110 = [s12 s23 s34 s49 s910 s101];
s39 = [s35 s56 s67 s78 s89];

T110 = [T12 T23 T34 T49 T910 T101];
T39 = [T35 T56 T67 T78 T89];
end