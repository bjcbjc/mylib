function RedGreenImage(data, mode, ylabels, xlabels)    %colormap(redgreencmap(256));        %image(colorscale(data, mode));    nanimagesc(colorscale(data, mode), genColorMap('rkg',255));    if nargin > 2        set(gca, 'YTick', 1:length(ylabels), 'YTickLabel',ylabels, 'FontSize',8);    end    if nargin > 3        set(gca, 'XTick', 1:length(xlabels), 'XTickLabel',xlabels, 'FontSize',8);    endendfunction [sc_data] = colorscale(data, mode)    switch mode        case 'overall'            M = max(max(data));  m = min(min(data));                    %sc_data = (256/(M-m))*(data-m);            sc_data = data;        case 'individual'            M = max(data');  m = min(data');            %sc_data = repmat(256./(M'-m'),1,size(data,2)) .* (data - repmat(m',1,size(data,2)));            globalM = max(M); globalm = min(m);                        sc_data = repmat((globalM-globalm)./(M'-m'),1,size(data,2))...                .* (data - repmat(m',1,size(data,2))) + globalm;        case '80percent'            %M = mean(prctile(data,90));            %m = mean(prctile(data,10));            M = prctile(reshape(data,1,numel(data)),90);            m = prctile(reshape(data,1,numel(data)),10);            MAbs = abs(M);            mAbs = abs(m);            if ((M>0 && m<0) && (MAbs > mAbs)) || (M>0 && m>0)                m=-M;            end            if ((M>0 && m<0) && (MAbs < mAbs)) || (M<0 && m<0)                M=-m;            end            %sc_data = (256/(M-m))*(data-m);            sc_data = data;            sc_data(data>M) = M;            sc_data(data<m) = m;        case 'rankind'            sc_data = zeros(size(data));            for i = 1:size(sc_data,1)                [tmp tmp sc_data(i,:)] = unique(data(i,:));                nNaN = sum(isnan(data(i,:)));                if nNaN > 0                    sc_data(i,end-(nNaN-1):end) = NaN;                end            end            M = max(sc_data'); m = min(sc_data');            globalM = max(M); globalm = min(m);            %sc_data = repmat(256./(M'-m'),1,size(sc_data,2)) .* (sc_data - repmat(m',1,size(sc_data,2)));            sc_data = repmat((globalM-globalm)./(M'-m'),1,size(sc_data,2))...                .* (sc_data - repmat(m',1,size(sc_data,2))) + globalm;        case 'rankall'            sc_data = zeros(size(data));            for i = 1:size(sc_data,1)                [tmp tmp sc_data(i,:)] = unique(data(i,:));                nNaN = sum(isnan(data(i,:)));                if nNaN > 0                    sc_data(i,end-(nNaN-1):end) = NaN;                end            end            %M = max(max(sc_data)); m = min(min(sc_data));            %sc_data = (256/(M-m))*(sc_data-m);    endend