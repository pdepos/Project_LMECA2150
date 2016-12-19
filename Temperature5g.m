function T5g = Temperature5g(lambda,h5g,x,y,z)

[mCO2_g,mH2O_g,mO2_g,mN2_g,Mg] = GasMassFraction(lambda,x,y,z);

h0g = mCO2_g*janaf('h','CO2',273.15) + mH2O_g*janaf('h','H2O',273.15)...
    + mO2_g*janaf('h','O2',273.15) + mN2_g*janaf('h','N2',273.15);

f = @(T)(mCO2_g*janaf('h','CO2',T) + mH2O_g*janaf('h','H2O',T)...
    + mO2_g*janaf('h','O2',T) + mN2_g*janaf('h','N2',T) - h0g - h5g);

opts = optimoptions(@fsolve,'Display','none');

T5g = fsolve(f,h5g+273.15,opts) - 273.15;

end