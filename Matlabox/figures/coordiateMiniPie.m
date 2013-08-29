function piehandles = coordiateMiniPie(X, Y, piedata, colors, piesize, pienames)
    % each column of piedata create a pie chart at location specified in X
    % and Y
    % piename: {index of pies, 'name'}
    %
    
    npie = size(piedata, 2);
    if npie ~= length(X) || npie ~= length(Y)
        error('dimension inconsistent');
    end
    if size(piedata,1) ~= size(colors,1)
        error('color dimension inconsistent');
    end
    
    if nargin < 5
        piesize = max( max(X)-min(X), max(Y)-min(Y) ) * 0.1;
    end
    if nargin < 6
        pienames = {[], []};
    end
    
    n = size(piedata,1);
    piestr = repmat({''}, 1, n);
    piehandles = cell(npie, 1);
    hold on
    for i = 1:npie
        piehandles{i} = pie(piedata(:, i), piestr);
        for j = 1:2:length(piehandles{i})
            set(piehandles{i}(j), 'vertices', ...
                bsxfun(@plus, get(piehandles{i}(j), 'vertices') .* piesize, [X(i), Y(i)]));
            set(piehandles{i}(j), 'facecolor', colors((j+1)/2, :), 'facealpha', 0.5, 'linestyle', 'none');
            set(piehandles{i}(j+1), 'position', ...
                get(piehandles{i}(j+1), 'position') .* piesize + [X(i), Y(i), 0]);
            if ismember(i, cell2mat(pienames(:,1))) && j == length(piehandles{i})-1
                [~, nameidx] = ismember(i, cell2mat(pienames(:,1)));
                set(piehandles{i}(j+1), 'string', pienames{nameidx, 2});
            end
        end        
    end
    hold off
    xl = xlim;
    yl = ylim;
    xlim([min(xl(1), yl(1)), max(xl(2), yl(2))]);
    ylim([min(xl(1), yl(1)), max(xl(2), yl(2))]);
end