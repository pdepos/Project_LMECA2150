function [s,h] = transitionHS(stateV,i,f,type,n)

if strcmp(type,'isobaric')
    
    s = linspace(stateV{i}.s,stateV{f}.s,n);
    h = zeros(1,n);
    h(1) = stateV{i}.h;
    h(n) = stateV{f}.h;
    for j = 2:n-1
        h(j) = XSteam('h_ps',stateV{i}.p,s(j));
    end
    
elseif strcmp(type,'saturation')
    
    s = linspace(stateV{i}.s,stateV{f}.s,n);
    h = zeros(1,n);
    h(1) = stateV{i}.h;
    h(n) = stateV{f}.h;
    for j = 2:n-1
        h(j) = XSteam('h_ps',stateV{i}.p,s(j));
    end
    
elseif strcmp(type,'undefined')
    
    s = linspace(stateV{i}.s,stateV{f}.s,n);
    p = linspace(stateV{i}.p,stateV{f}.p,n);
    h = zeros(1,n);
    h(1) = stateV{i}.h;
    h(n) = stateV{f}.h;
    for j = 2:n-1
        h(j) = XSteam('h_ps',p(j),s(j));
    end
    
else
    s = zeros(n,1);
    h = zeros(n,1);
    disp('(h,s) diagram error : unknown type of transition');
end

end