function [boxplotHandle, dataHandle] = boxplotWithData(data, group, varargin)
para.marker = '.';
para.markersize = 14;
para.markercolor = ColorDict.translate('n');
para.jitter = 0.3;
para.withinGroupLabel = '';

passthrough = varargin;
[~, i] = ismember(fieldnames(para), varargin(1:2:end));
i(i==0) = [];
if ~isempty(i)
    i = i(:)*2-1;
    paraidx = reshape([i,i+1]', length(i)*2, 1);
    paradata = varargin( paraidx );
    para = assignpara(para, paradata{:});
    passthrough = varargin(setdiff(1:length(varargin), paraidx));
end

if isempty(para.withinGroupLabel)
    para.withinGroupLabel = ones(size(data));
else
    if size(para.withinGroupLabel, 2) == 1 && size(data,2) > 1
        para.withinGroupLabel = repmat(para.withinGroupLabel, 1, size(data,2));
    end
end

horizontal = ismember('horizontal', passthrough(cellfun(@ischar,passthrough)));
if nargin < 2
    group = repmat(1:size(data,2), size(data, 1), 1);
end

if any(size(group) ~= size(data))
    group = group(:)';
    if size(group, 2) == length(group)        
        group = repmat(group, size(data,1), 1);
    else
        error('cannot proceed');
    end
end

if size(data,2)>1
    data = data(:);
    group = group(:);
    para.withinGroupLabel = para.withinGroupLabel(:);
end

if isnumeric(group)
    [ugroup, ~, groupIdx] = unique(group);
else
    [ugroup, ~, groupIdx] = unique(group, 'stable');
end


[uWithinGroupLabel, ~, withinGroupLabelIdx] = unique(para.withinGroupLabel);
nWithinGroupLabel = length(uWithinGroupLabel);
if nWithinGroupLabel > 1 && size(para.markercolor,1) == 1
    para.markercolor = ColorDict.translate({'n', 'p', 'g', 'o', 'c'});
end
h = NaN(length(ugroup),nWithinGroupLabel);
cla
hold all;
if nWithinGroupLabel == 1    
    for i = 1:length(ugroup)    
        idx = find(groupIdx == i);
        if horizontal
            h(i) = plot(data(idx), i + (rand(length(idx),1)-0.5)*para.jitter, ...
                para.marker, 'color', para.markercolor, 'markersize', para.markersize);
        else
            h(i) = plot(i + (rand(length(idx),1)-0.5)*para.jitter, data(idx), ...
                para.marker, 'color', para.markercolor, 'markersize', para.markersize);
        end    
    end
else    
    nColor = size(para.markercolor, 1);
    for i = 1:length(ugroup)    
%         idx = find(groupIdx == i);
        for j = 1:nWithinGroupLabel
            idx = find(groupIdx == i & withinGroupLabelIdx == j);
            if horizontal
                h(i,j) = plot(data(idx), i + (rand(length(idx),1)-0.5)*para.jitter, ...
                    para.marker, 'color', para.markercolor(mod(j+1,nColor)+1,:), 'markersize', para.markersize);
            else
                h(i,j) = plot(i + (rand(length(idx),1)-0.5)*para.jitter, data(idx), ...
                    para.marker, 'color', para.markercolor(mod(j+1,nColor)+1,:), 'markersize', para.markersize);
            end    
        end
    end
    for j = 1:nWithinGroupLabel
        set(h(1,j), 'displayname', uWithinGroupLabel{j});
    end
    legend(h(1,:));
end

if ~ismember('symbol', passthrough(1:2:end))
    passthrough(end+1:end+2) = {'symbol', ''};
end
bh = boxplot(data, group, passthrough{:});
hold off;

% set(gca, 'ytick', 1:length(ugroup), 'yticklabel', ugroup);
if horizontal
    set(gca, 'xtickmode', 'auto', 'xticklabelmode', 'auto');
else
    set(gca, 'ytickmode', 'auto', 'yticklabelmode', 'auto');
end
if nargout > 0
    boxplotHandle = bh;
end
if nargout > 1
    dataHandle = h;
end
