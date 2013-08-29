function rankdata = rankData(data, mode, tiekeep, nakeep)
    %rank samples in each series of data
    %data: #sample x #series
    %mode: 'descend' or 'ascend' (default)
    %tiekeep: true if keep tie; otherwise, ties will be ranked randomly
    %nakeep: true if keep NaN; otherwise, NaN will be ranked as large
    %numbers
    %    
    %
    
    if nargin < 2, mode = 'ascend'; end
    if nargin < 3, tiekeep = false; end
    if nargin < 4, nakeep = false; end
    
    [nsample, nseries] = size(data);    
    rankdata = zeros(nsample, nseries);
    
    if ~tiekeep && ~nakeep
        for i = 1:nseries
            [~, si] = sort(data(:,i),mode);
            rankdata(si,i) = (1:nsample)';
        end
    elseif ~tiekeep && nakeep
        for i = 1:nseries
            nani = isnan(data(:,i));
            [~, si] = sort(data(~nani,i),mode);
            j = find(nani==0);
            rankdata(j(si),i) = (1:length(j))';
            rankdata(nani,i) = NaN;
        end
    elseif tiekeep && ~nakeep
        tmp = data;
        tmp(isnan(tmp)) = Inf; %NaN is not unique
        sirank = NaN(nsample, 1);
        for i = 1:nseries
            [~, si] = sort(tmp(:,i), mode);            
            sirank(si) = 1:nsample;
            u = unique(tmp(:,i));
            for j = 1:length(u)
                v = tmp(:,i) == u(j);
                rankdata(v, i) = mean(sirank(v));
            end
        end
%         for i = 1:nseries
%             u = sort(unique(tmp(:,i)), mode);            
%             for j = 1:length(u)
%                 rankdata(tmp(:,i)==u(j),i) = j;
%             end
%         end
    else %tiekep and nakeep
        sirank = NaN(nsample, 1);
        for i = 1:nseries
            vi = find(~isnan(data(:,i)));
            [~, si] = sort(data(vi, i), mode);            
            sirank(vi(si)) = 1:length(vi);
            u = unique(data(vi, i));
            for j = 1:length(u)
                v = data(:,i) == u(j);
                rankdata(v, i) = mean(sirank(v));
            end
        end
%         for i = 1:nseries
%             nani = isnan(data(:,i));
%             u = sort(unique(data(~nani,i)),mode);
%             for j = 1:length(u)
%                 rankdata(data(:,i)==u(j),i) = j;
%             end
%             rankdata(nani,i) = NaN;
%         end
    end
end