function [counts, uniqElements, uindex] = eleCounts(data, sortByCount)%, group)
    % a basic function that counts the occurence of elements in a vector
    % return unique elements in uniqElements
    % return counts corresponding to elements in uniqElements
    %
    % data: vector or cell; if data is a vector, each element is a number; 
    %   if data is cell array, an element is a cell.    
    % group: a vector indicating the group of each vector; if omitted, each
    %   element in data is assumed to be a different group
    %
    % If multiple groups are specified, the counting is done to see how
    % many "different groups" the element appears. If the element appears
    % more than once in the same group, it will only be counted once.
    %
    % return:
    %   counts: counts for each unique element
    %   uniqElements: unique elements
    %   grmembers: groups that belong to each unique element
    %
    
    if nargin < 2, sortByCount = true; end
    if nargout == 3, sortByCount = false; end
    data = data(:);
    [uniqElements, ~, uindex] = unique(data);
    n = length(uniqElements);
    counts = zeros(n,1);
    for i = 1:n
        counts(i) = sum(uindex==i);
    end
    if sortByCount
        [~,si] = sort(counts, 'descend');
        counts = counts(si);
        uniqElements = uniqElements(si);
    end
    
%     data = data(:);
% 
%     if nargin < 2
%         group = 1:length(data);
%     else
%         if length(group) ~= length(data)
%             error('number of group ~= number of data\n');
%         end
%     end
%     
%             
%     if iscell(data)        
%         if iscellstr(data)
%             extvec = data;
%             extgroup = group;
%         else
%             extvec = zeros(0,1);
%             extgroup = zeros(0,1);
%             unigroup = unique(group);
%             for i = 1:length(unigroup)
%                 tmp = [];
%                 subdata = data(group==unigroup(i));
%                 for j = 1:length(subdata)
%                     tmp = union(tmp, subdata{j});
%                 end
%                 [nrow ncol] = size(tmp);
%                 if nrow < ncol
%                     tmp = tmp';
%                 end
%                 extvec = [extvec; tmp];            
%                 extgroup = [extgroup; repmat(unigroup(i),length(tmp),1)];
%             end
%         end
%     else
%         extvec = data;        
%         extgroup = group;
%     end
%     
%     
%     if ~iscellstr(data)
%         uniqElements = unique(extvec);
%         counts = zeros(size(uniqElements));
%         grmembers = cell(size(uniqElements));
%         for i = 1:length(uniqElements)
%             mask = extvec == uniqElements(i);
%             counts(i) = length(unique(extgroup(mask)));
%             grmembers{i} = unique(extgroup(mask));
%         end
%     else
%         uniqElements = unique(extvec);
%         counts = zeros(size(uniqElements));
%         grmembers = cell(size(uniqElements));
%         for i = 1:length(uniqElements)
%             mask = strmatch(uniqElements{i}, extvec, 'exact');
%             counts(i) = length(unique(extgroup(mask)));
%             grmembers{i} = unique(extgroup(mask));
%         end
%     end
% 
%     %sort by counts
%     [tmp si] = sort(counts, 'descend');
%     counts = counts(si);
%     uniqElements = uniqElements(si);
%     grmembers = grmembers(si);
end