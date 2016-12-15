% Script calculating the points of the oswald diagram and the different
% flame temperatures. 

qair   = 18.13;   % Modify this variable to study a different situation.


%% Oswald Diagram

P = 0.11914;
Q = 0.21;
S = 0.056214;
alpha1 = -P/S;
alpha2 = -P/Q;
beta = P;
x = linspace(0,0.25);
y1 = alpha1*x + beta;
y2 = alpha2*x + beta;

L  = [20/18.13 23/18.13 26/18.13];

CO2  = [1 1 1];
O2   = [1 1 1];
for i = 1:3
    CO2(i) = 1.0156 / (1.0156 + 1.96435*(L(i) - 1) + (7.386*L(i) + 0.123));
    O2(i)  = 1.96435*(L(i) - 1) / (1.0156 + 1.96435*(L(i) -1) + (7.386*L(i)...
    + 0.123));
end

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(x,y1,x,y2,O2(1),CO2(1),'*',O2(2),CO2(2),'*',O2(3),CO2(3),'*');
set(plot1(1),'DisplayName','? = 1','Color',[1 0 0]);
set(plot1(2),'DisplayName','k = 0','Color',[0 0 1]);

% Create xlabel
xlabel({'[O_2]'''},'FontSize',11);

% Create title
title({'Oswald Diagram'},'FontSize',11);

% Create ylabel
ylabel({'[CO_2]'''},'FontSize',11);

% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0 0.13]);
box(axes1,'on');
% Create legend
legend(axes1,'show');

% Create textbox
annotation(figure1,'textbox',...
    [0.130483689538808 0.833514769058355 0.0583657575952404 0.0599999988741345],...
    'String',{'P'},...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.308211473565804 0.0979310841075814 0.0583657575952404 0.0599999988741345],...
    'String',{'S'},...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.778402699662544 0.100744023629382 0.0622568080629357 0.0599999988741345],...
    'String',{'Q'},...
    'EdgeColor','none');

%% Flame Temperature Calculations
MWO2   = 31.998;
MWN2   = 28.014;
MWCO2  = 44.008;
MWH2O  = 18.01494;
MWH2   = 2.01594;
MWC    = MWCO2 - MWO2;
MWCH4  = MWC + 2*MWH2;
MWCH18 = MWC + 1.8 / 2 * MWH2;
qgas   = 2.05; %m^3/h
qair0  = 18.13;
MWgas  = 18.6244;
MWair  = 0.21*MWO2 + 0.79*MWN2;
VM     = 24.3; %m^3 / kmol
lambda = qair/qair0;
Ta     = 23;
mair   = qair*MWair/VM;
mgas   = qgas*MWgas/VM;
mfumes = mair + mgas;

Wfgasses  = (1.0156*MWCO2 + 1.9287*MWH2O + (lambda - 1)*1.96435*MWO2 + ...
            (7.386*lambda + 0.123)*MWN2);             

CO2_frac = 1.0156*MWCO2/Wfgasses;        %kg_CO2 / kg_fg
H2O_frac = 1.9287*MWH2O/Wfgasses;        %kg_H2O / kg_fg
O2_frac  = (lambda - 1)*1.96435*MWO2/Wfgasses;
N2_frac  = (7.386*lambda + 0.123)*MWN2/Wfgasses;

tfg0 = 1500;
cpO200  = janaf('c','O2',273.15 + Ta);
cpO2f0  = janaf('c','O2',tfg0);
cpN200  = janaf('c','N2',273.15 + Ta);
cpN2f0  = janaf('c','N2',tfg0);
cpCO200 = janaf('c','CO2',273.15 + Ta);
cpCO2f0 = janaf('c','CO2',tfg0);
cpH2O00 = janaf('c','H2O',273.15 + Ta);
cpH2Of0 = janaf('c','H2O',tfg0);

% cp for flue gasses first guess
cpfg0 = 0.5 * (CO2_frac*(cpCO200 + cpCO2f0) + H2O_frac*(cpH2O00 + cpH2Of0) ...
             + O2_frac*(cpO200  +  cpO2f0) +  N2_frac*(cpN200  +  cpN2f0)); %kJ / kgfg K

         
tfgprev = tfg0;
tfgnew  = 293 + 63000/ (mfumes*cpfg0);
while abs(tfgprev - tfgnew) > 2
    cpO20  = janaf('c','O2',273.15 + Ta);
    cpO2f  = janaf('c','O2',tfgnew);
    cpN20  = janaf('c','N2',273.15 + Ta);
    cpN2f  = janaf('c','N2', tfgnew);
    cpCO20 = janaf('c','CO2',273.15 + Ta);
    cpCO2f = janaf('c','CO2',tfgnew);
    cpH2O0 = janaf('c','H2O',273.15 + Ta);
    cpH2Of = janaf('c','H2O',tfgnew);

    % cp for flue gasses
    cpfgnew = 0.5 * (CO2_frac*(cpCO20 + cpCO2f) + H2O_frac*(cpH2O0 + cpH2Of) ...
             + O2_frac*(cpO20  +  cpO2f) +  N2_frac*(cpN20  +  cpN2f)); %kJ / kgfg K
    tfgprev = tfgnew;
    tfgnew  = 293 + 63000 / (mfumes*cpfgnew);
end
Tad = tfgnew - 273.15;
cpf = cpfgnew;


%% Measured Temperature Profiles
xd1   = [5 5 5 5 5 5 5];
xd2   = [10 10 10 10 10 10 10 10 10 10 10 10 10];
xd3   = [15 15 15 15 15 15 15 15 15 15 15 15 15];
yd1   = [-3 -2 -1 0 1 2 3];
yd    = [-6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6];
tc120 = [357 514 712 765 712 514 357];
tc220 = [566 600 652 670 760 783 792 783 760 670 652 600 566];
tc320 = [685 723 782 782 818 831 845 831 818 782 782 723 685];
tc123 = [353 460 705 811 705 460 353];
tc223 = [548 573 615 654 695 697 702 697 695 654 615 573 548];
tc323 = [663 688 730 768 776 787 790 787 776 768 730 688 663];
tc126 = [370 465 690 785 690 465 370];
tc226 = [510 547 581 606 652 680 672 680 652 606 581 547 510];
tc326 = [605 628 670 728 750 764 760 764 750 728 670 628 605];

% Uncomment to get 2D plots of the profiles
% figure
% subplot(2,1,1)
% plot(yd,tc220,'or',yd,tc320,'*r',...
%     yd,tc223,'ob',yd,tc323,'*b',yd,tc226,'og',yd,tc326,'*g');
% subplot(2,1,2)
% plot(yd1,tc120,'xr',yd1,tc123,'xb',yd1,tc126,'xg');

% 3D plot of the temperature profiles
plot3(xd1,yd1,tc120,'r',xd1,yd1,tc123,'b',xd1,yd1,tc126,'g',xd2,yd,tc220,'r',...
    xd2,yd,tc223,'b',xd2,yd,tc226,'g',xd3,yd,tc320,'r',xd3,yd,tc323,'b',xd3,...
    yd,tc326,'g');
