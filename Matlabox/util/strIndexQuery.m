function [indices] = strIndexQuery(cells, queries)
%this function only return the first one searched
%that means if more than one entry is in the data pool,
%it will only return the first occurence

    if ~iscell(queries)
        queries = {queries};
    end
    %t = cputime();
    indices = [];
    for i = 1:length(queries)
%         found = 0;
%         for j = 1:length(cells)
%             if strcmp(cells{j},queries{i})
%                 indices = [indices, j];
%                 found = 1;
%                 break
%             end
%         end
%         if found == 0
%             indices = [indices, -1];
%         end
        j = strmatch(queries{i}, cells, 'exact');
        if isempty(j)
            indices = [indices; -1];
        else
            indices = [indices; j(1)];
        end
    end
    %fprintf('time:%f\n',cputime()-t);
end