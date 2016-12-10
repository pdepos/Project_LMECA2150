function T5g = Temperature5g(lambda,h5g,x,y,z)

[mCO2_g,mH2O_g,mO2_g,mN2_g,Mg] = GasMassFraction(lambda,x,y,z);

h0g = mCO2_g*janaf('h','CO2',273.15) + mH2O_g*janaf('h','H2O',273.15)...
    + mO2_g*janaf('h','O2',273.15) + mN2_g*janaf('h','N2',273.15);

f = @(T)(mCO2_g*janaf('h','CO2',T) + mH2O_g*janaf('h','H2O',T)...
    + mO2_g*janaf('h','O2',T) + mN2_g*janaf('h','N2',T) - h0g - h5g);

%%%% initial points %%%%
a = 350;
if f(a)>=0
    while f(a) >=0
       a = a-10; 
    end
end

b = 650;
if f(b)<=0
   while f(b)<=0
       b = b+10;
   end
end

%%%% iterative solver %%%%
while (b-a)>0.0001
    c = (a+b)/2;
    if sign(f(c)) ~= sign(f(a))
        b = c;
    else
        a = c;
    end
end

T5g = ((a+b)/2) - 273.15; 
end