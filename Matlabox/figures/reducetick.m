function reducetick(axhandle, maxntick, xory)
    % maxntick: max number of ticks wanted
    % xory: 1x2, specify reduce x and/or y axis
    tickname = {'xtick', 'ytick'};    
    maxntick = maxntick - 2; %keep the first and last ticks
    reducerange = 2:maxntick;
    for i = 1:length(xory)
        newticks = [];
        if xory(i) ~= 0
            ticks = get(axhandle, tickname{i});
            curntick = length(ticks);
            if curntick > (maxntick+2)
                for j = 1:length(reducerange)
                    if mod(curntick-1, reducerange(j)) == 0
                        newticks = ticks(1:(curntick-1)/reducerange(j):end);
                        break
                    end
                end
                if ~isempty(newticks)
                    set(axhandle, tickname{i}, newticks);
                end
            end
        end
    end