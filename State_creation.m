function State = State_creation(n,alpha,beta)
% function creating an cell vecto 'State' of n lignes where n is 
% the number of cycle's states.
% Each cell contain a structure with the following fields :
%   - T : the temperature in °C
%   - p : the pressure in bar
%   - h : the enthalpy in kJ/kg
%   - s : the entropy in kJ/(kg*°C)
%   - x : the title (no unit)
% The differents fields have a default value of 0 except the title x for
% which the default value is NaN

State = cell(n,1);

for i = 1:n
    if i <= alpha
        formatSpec = '%d';
        str = sprintf(formatSpec,i);
        State{i} = struct('States',str,'t',0,'p',0,'x',NaN,'h',0,'s',0,'e',0);
    else 
        formatSpec = '%d.%d';
        modulo = mod(i,beta);
        if modulo == 0
            a = 4;
        elseif modulo == 1
            a = 6;
        elseif modulo == 2
            a = 7;
        end
        b = floor(i./beta) - 2;
        str = sprintf(formatSpec,a,b);
        State{i} = struct('States',str,'t',0,'p',0,'x',NaN,'h',0,'s',0,'e',0);
    end
end
end