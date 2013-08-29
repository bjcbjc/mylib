function dist = distRegion(region1, region2, mode)
    %calculate the distance between two regions
    %region1, region2: marker strings or matrix of two columns (start,
    %   end; start <= end)
    %mode: {'min'}
    %
    
    if nargin < 3,  mode = 'min';   end
    
    region1 = parseLoc(region1);
    region2 = parseLoc(region2);
    nloc1 = size(region1,1);
    nloc2 = size(region2,1);
    
    dist = NaN(nloc1, nloc2);
    if strcmp(mode, 'min')
        meanloc1 = mean(region1,2);
        meanloc2 = mean(region2,2);
        for i = 1:nloc1
            j = meanloc2 <= meanloc1(i); %those in reg2 are in front of region1(i)
            dist(i,j) = min(region1(i,:)) - max(region2(j,:),[],2);
            dist(i,~j) = min(region2(~j,:),[],2) - max(region1(i,:));            
        end
        dist(dist<0) = 0;
    end
    
end

function loc = parseLoc(locstr)
    if ischar(locstr) %single marker string
        [chrm ps pe] = parseMarker(locstr);
        loc = [min([ps pe]) max([ps pe])];
    elseif iscellstr(locstr) %multiple marker strings
        nm = length(locstr);
        loc = NaN(nm,2);
        for i = 1:nm
            [chrm ps pe] = parseMarker(locstr{i});
            loc(i,1) = min([ps pe]);
            loc(i,2) = max([ps pe]);
        end
    elseif isnumeric(locstr) && size(locstr,2) == 2
        loc = locstr;
    else
        error('Unknown input.\n');
    end
end

function [chrm ps pe] = parseMarker(marker)
    [chrm remain] = strtok(marker,'_');
    chrm = str2num(strrep(chrm,'M',''));
    [ps pe] = strtok(remain,'_');
    ps = str2num(ps);
    pe = str2num(strrep(pe,'_',''));
end