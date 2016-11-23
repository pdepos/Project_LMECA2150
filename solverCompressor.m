function [T2,na] = solverCompressor(T2_0,etaPiC,Ra,T1,r,mO2_a,mN2_a)

% step 1

if T1 < 300
    n = round(300-T1);
    T = [ones(1,n)*300 300:1:T2_0];
    cp = mean(mO2_a*janaf('c','O2',T) + mN2_a*janaf('c','N2',T))*10^3;
    na = (1 - (Ra/(etaPiC*cp)))^(-1); 
    T2 = T1*(r^((na-1)/na));
else
    T = [T1:1:T2_0];
    cp = mean(mO2_a*janaf('c','O2',T) + mN2_a*janaf('c','N2',T))*10^3;
    na = (1 - (Ra/(etaPiC*cp)))^(-1);
    T2 = T1*(r^((na-1)/na));
end

T2_prev = T2_0;
ma_prev = 5;

while (abs(T2-T2_prev) > 0.0001) || (abs(na-ma_prev) > 0.0001)
    T2_prev = T2;
    ma_prev = na;
    if T1 < 300
        n = round(300-T1);
        T = [ones(1,n)*300 300:1:T2];
        cp = mean(mO2_a*janaf('c','O2',T) + mN2_a*janaf('c','N2',T))*10^3;
        na = (1 - (Ra/(etaPiC*cp)))^(-1);
        T2 = T1*(r^((na-1)/na));
    else
        T = [T1:1:T2];
        cp = mean(mO2_a*janaf('c','O2',T) + mN2_a*janaf('c','N2',T))*10^3;
        na = (1 - (Ra/(etaPiC*cp)))^(-1);
        T2 = T1*(r^((na-1)/na));
    end
end


end