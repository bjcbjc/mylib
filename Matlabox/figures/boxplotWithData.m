function [boxplotHandle, dataHandle] = boxplotWithData(data, group, varargin)
para.marker = '.';
para.markersize = 14;
para.markercolor = ColorDict.translate('n');
para.jitter = 0.3;

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

horizontal = ismember('horizontal', passthrough(cellfun(@ischar,passthrough)));
if size(data,2)>1 && nargin < 2
    group = repmat(1:size(data,2), size(data, 1), 1);
    data = data(:);
    group = group(:);
end

if isnumeric(group)
    [ugroup, ~, groupIdx] = unique(group);
else
    [ugroup, ~, groupIdx] = unique(group, 'stable');
end
h = NaN(length(ugroup),1);
cla
hold on;
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
