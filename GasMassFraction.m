function [mCO2_g,mH2O_g,mO2_g,mN2_g,Mg] = GasMassFraction(lambda,x,y,z)

% molar fractions
moltot = x + y/2 + (lambda-1)*(x+(y-2*z)/4) + 3.76*lambda*(x+(y-2*z)/4);
molCO2_g = x/moltot;
molH2O_g = (y/2)/moltot;
molO2_g = ((lambda-1)*(x+(y-2*z)/4))/moltot;
molN2_g = (3.76*lambda*(x+(y-2*z)/4))/moltot;

% molar mass
Mg = molCO2_g*44 + molH2O_g*18 + molO2_g*32 + molN2_g*28;

% mass fractions
mCO2_g = molCO2_g*(44/Mg);
mH2O_g = molH2O_g*(18/Mg);
mO2_g = molO2_g*(32/Mg);
mN2_g = molN2_g*(28/Mg);


end
    
