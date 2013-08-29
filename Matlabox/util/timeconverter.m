function res = timeconverter(timestr, unit)
    %convert time string to second (or minute)
   
    if nargin < 2, unit = 'min'; end
    
    if ischar(timestr)
        timestr = {timestr};
    end
    
    n = length(timestr);
    

    
    if strcmpi(unit, 'hour') || strcmpi(unit, 'hr')
        uplace = 1;
    elseif strcmpi(unit, 'min') || strcmpi(unit, 'm')
        uplace = 2;
    else
        uplace = 3;
    end
        
    h = zeros(n, 3);
    for i = 1:n
        tformat = length(strfind(timestr{i}, ':'));
        if tformat == 1
            h(i,2:3) = sscanf(timestr{i}, '%d:%d'); %hour and min
        elseif tformat == 2
            h(i,:) = sscanf(timestr{i}, '%d:%d:%d'); %h:m:s
        else
            error('unknown format\n');
        end
    end
    
    res = zeros(n, 1);
    for i = 1:3
        res = res + h(:,i)*(60^(uplace-i));
    end
end