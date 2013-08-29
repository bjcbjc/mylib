function newstr = breaklongstr(s, len, breaker, onlyafter)

    if nargin < 4, onlyafter = ''; end

    ntest = length(onlyafter);
    
    newstr = '';    
    addbreaker = false;
    for i = 1:length(s)
        if mod(i, len) == 0
            addbreaker = true;            
        end
        if addbreaker
            if isempty(onlyafter)                
                %newstr = strcat(newstr, breaker);
                newstr = sprintf('%s%s', newstr, breaker);
                addbreaker = false;
            elseif strcmpi(newstr(end-ntest+1:end), onlyafter)                
                %newstr = strcat(newstr, breaker);
                newstr = sprintf('%s%s', newstr, breaker);
                addbreaker = false;
            end
        end
        newstr = sprintf('%s%s', newstr, s(i));
    end