function [s,T] = transitionTS(stateV,i,f,type,n)

if strcmp(type,'isobaric')
    
    s = linspace(stateV{i}.s,stateV{f}.s,n);
    T = zeros(1,n);
    T(1) = stateV{i}.T;
    T(n) = stateV{f}.T;
    for j = 2:n-1
        T(j) = XSteam('T_ps',stateV{i}.p,s(j));
    end
    
elseif strcmp(type,'saturation')
    
    s = linspace(stateV{i}.s,stateV{f}.s,n);
    T = ones(1,n)*stateV{i}.T;
    
elseif strcmp(type,'undefined')
    
    h = linspace(stateV{i}.h,stateV{f}.h,n);
    s = linspace(stateV{i}.s,stateV{f}.s,n);
    T = zeros(1,n);
    T(1) = stateV{i}.T;
    T(n) = stateV{f}.T;
    for j = 2:n-1
        T(j) = XSteam('T_hs',h(j),s(j));
    end
    
else
    s = zeros(n,1);
    T = zeros(n,1);
    disp('(T,s) diagram error : unknown type of transition');
end



end