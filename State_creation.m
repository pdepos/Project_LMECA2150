function State = State_creation(n,alpha,beta)
% function creating an cell vecto 'State' of n lignes where n is 
% the number of cycle's states.
% Each cell contain a structure with the following fields :
%   - T : the temperature in �C
%   - p : the pressure in bar
%   - h : the enthalpy in kJ/kg
%   - s : the entropy in kJ/(kg*�C)
%   - x : the title (no unit)
% The differents fields have a default value of 0 except the title x for
% which the default value is NaN

State = cell(n,1);
for i = 1:n
    if i <= 2
        formatSpec1 = '%d';
        str1 = sprintf(formatSpec1,i);
        State{i} = struct('States',str1,'t',0,'p',0,'x',NaN,'h',0,'s',0,'e',0);
    elseif i == 3
        formatSpec1 = '%d''';
        str1 = sprintf(formatSpec1,2);
        State{i} = struct('States',str1,'t',0,'p',0,'x',NaN,'h',0,'s',0,'e',0);
    elseif i == 4
        formatSpec1 = '%d%s';
        str1 = sprintf(formatSpec1,2,'"');
        State{i} = struct('States',str1,'t',0,'p',0,'x',NaN,'h',0,'s',0,'e',0);
    elseif (i > 4) && ( i <= alpha)
        formatSpec1 = '%d';
        str1 = sprintf(formatSpec1,i-2);
        State{i} = struct('States',str1,'t',0,'p',0,'x',NaN,'h',0,'s',0,'e',0);
    elseif i > alpha
        formatSpec1 = '%d.%d';
        modulo = mod((i-2),beta);
        if modulo == 1
            a = 4;
        elseif modulo == 2
            a = 6;
        elseif modulo == 3
            a = 7;
        elseif modulo == 0
            a = 8;
        end
        roman_index = floor((i-3)./beta)-1;
        str1 = sprintf(formatSpec1,a,roman_index);
        State{i} = struct('States',str1,'t',0,'p',0,'x',NaN,'h',0,'s',0,'e',0);
    end
end
end