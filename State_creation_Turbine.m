function State = State_creation_Turbine(n)
% function creating an cell vecto 'State' of n lignes where n is 
% the number of cycle's states.
% Each cell contain a structure with the following fields :
%   - T : the temperature in �C
%   - p : the pressure in bar
%   - h : the enthalpy in kJ/kg
%   - s : the entropy in kJ/(kg*�C)
% The differents fields have a default value of 0 except the title x for
% which the default value is NaN

State = cell(n,1);

for i = 1:n
    State{i} = struct('T',0,'p',0,'h',0,'s',0,'e',0);
end
end