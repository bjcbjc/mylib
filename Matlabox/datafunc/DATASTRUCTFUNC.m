classdef DATASTRUCTFUNC < handle
    
    methods (Static = true)
        % query data.(keyfield) for "queries" and return all data in
        % queryFields with matching queries
        function [result, idx] = queryData(data, keyfield, queries, queryFields)
            
            if ischar(queryFields)
                queryFields = {queryFields};
            end
            if ~ismember(keyfield, fieldnames(data))
                fprintf('key field is not in data\n');
                result = [];
                return
            end
            if any(~ismember(queryFields, fieldnames(data))) 
                fprintf('some fields are not in data\n');
                result = [];
                return
            end
            nfields = length(queryFields);
            [~, idx] = ismember(queries, data.(keyfield));
            
            n = length(queries);
            keyDim = length(data.(keyfield));
            if nfields == 1
                [dim1, dim2] = size(data.(queryFields{1}));
                %which dim matches key
                if keyDim == dim1
                    result = DATASTRUCTFUNC.allocateData(class(data.(queryFields{1})), [n, dim2]);
                    result(idx~=0, :) = data.(queryFields{1})(idx(idx~=0 ), :);
                elseif keyDim == dim2
                    result = DATASTRUCTFUNC.allocateData(class(data.(queryFields{1})), [dim1, n]);
                    result(:, idx~=0) = data.(queryFields{1})(:, idx(idx~=0 ));
                else
                    error('cannot determine which dimension to match for %s', queryFields{1});
                end                                
            else
                result = cell(n, nfields);
                result(:) = {'NA'};
                allnumeric = true;
                for i = 1:length(queryFields)                    
                    [dim1, dim2] = size(data.(queryFields{i}));
                    if keyDim == dim1
                        if dim2 > 1
                            fprintf('%s is a matrix, take only first column\n', queryFields{i});
                        end
                        dataForCopy = data.(queryFields{i})(idx(idx ~=0 ), 1);                        
                    elseif keyDim == dim2
                        if dim1 > 1
                            fprintf('%s is a matrix, take only first row\n', queryFields{i});
                        end
                         dataForCopy = data.(queryFields{i})(1, idx(idx ~= 0));                        
                    else
                        error('cannot determine which dimension to match for %s', queryFields{i});
                    end
                    if iscell(data.(queryFields{i}))
                        allnumeric = false;
                        if keyDim == dim1
                            result(idx~=0, i) = dataForCopy;
                        else
                            result(i, idx~=0) = dataForCopy;
                        end
                    else
                        if keyDim == dim1
                            result(idx~=0, i) = num2cell(dataForCopy);
                        else
                            result(i, idx~=0) = num2cell(dataForCopy);
                        end
                    end                    
                end
                if allnumeric
                    result = cell2mat(result);
                end
            end
        end
        
        function data = copyData(data, copyFrom, dataKey, copyKey, copyFields, newFieldNames)
            if nargin < 6, newFieldNames = {}; end
            fds = fieldnames(copyFrom);
            existfds = fieldnames(data);
            if ~iscell(copyFields)
                copyFields = {copyFields};                
            end
            if ~iscell(newFieldNames)
                newFieldNames = {newFieldNames};
            end
            if isempty(newFieldNames)
                newFieldNames = copyFields;
            else
                if length(newFieldNames) ~= length(copyFields)
                    error('#newFieldNames ~= #copyFields');
                end
            end
            [~, idx] = ismember(data.(dataKey), copyFrom.(copyKey));
            n = length(idx);
            for i = 1:length(copyFields)
                if ~ismember(copyFields{i}, fds)
                    fprintf('skip %s, it is not in data\n', copyFields{i});
                    continue
                end
                if ismember(newFieldNames{i}, existfds)
                    fprintf('%s exists, create %s_2\n', newFieldNames{i}, newFieldNames{i});
                    newfd = sprintf('%s_2', newFieldNames{i});
                else
                    newfd = newFieldNames{i};
                end
                [dim1, dim2] = size(copyFrom.(copyFields{i}));
                %which dim matches key
                if length(copyFrom.(copyKey)) == dim1
                    data.(newfd) = DATASTRUCTFUNC.allocateData(class(copyFrom.(copyFields{i})), [n, dim2]);
                    data.(newfd)(idx ~= 0, :) = copyFrom.(copyFields{i})( idx(idx~=0), : );
                elseif length(copyFrom.(copyKey)) == dim2
                    data.(newfd) = DATASTRUCTFUNC.allocateData(class(copyFrom.(copyFields{i})), [dim1, n]);
                    data.(newfd)(:, idx ~= 0) = copyFrom.(copyFields{i})( :, idx(idx~=0));
                else
                    error('cannot determine which dimension to match for %s', copyFields{i});
                end
            end
        end
            
        % order all fields in datastruct according to datastruct.(keyfield)
        % and orderBy (new order)
        function newdb = orderDataStruct(datastruct, keyfield, orderBy)
            % input:
            %   dbstruct, structure containing the database
            %   fd2align: string, the fieldname to align
            %   alignlist: cell array, list to align to
            %
            % Note: for those that are not in alignlist, create NaN in the matrix
            %
            [~, i] = ismember(orderBy, datastruct.(keyfield));
            
            if all(i==0)
                error('nothing matched! Use your brain!');
            end
            
            fds = fieldnames(datastruct);
            nfd = length(fds);
            dim2change = zeros(nfd,1);
            ndim = zeros(nfd,1);
            
            matchdim = length(datastruct.(keyfield));
            for j = 1:nfd
                dim = size(datastruct.(fds{j}));
                ndim(j) = length(dim);
                for di = 1:ndim(j)
                    if dim(di) == matchdim
                        dim2change(j) = di;
                        break
                    end
                end
            end
            
            newlen = length(orderBy);
            for j = 1:nfd
                if dim2change(j) == 0
                    newdb.(fds{j}) = datastruct.(fds{j});
                else
                    if strcmp(fds{j}, keyfield)
                        newdb.(fds{j}) = orderBy;
                    else
                        ifsparse = issparse(datastruct.(fds{j}));
                        permorder = [dim2change(j) setdiff(1:ndim(j), dim2change(j))];
                        backorder = NaN(1, ndim(j));
                        backorder(dim2change(j)) = 1;
                        backorder(setdiff(1:ndim(j), dim2change(j))) = 2:ndim(j);
                        d = permute(datastruct.(fds{j}), permorder);
                        dim = size(datastruct.(fds{j}));
                        newdim = dim(permorder);
                        newdim(1) = newlen;
                        if iscell(d)
                            newdb.(fds{j}) = cell(newdim);
                        else
                            newdb.(fds{j}) = NaN(newdim);
                        end
                        s = '';
                        for k = 1:ndim-1
                            s = [s ',:'];
                        end
                        s = [s ')'];
                        eval(sprintf('newdb.(fds{j})(i~=0%s = d(i(i~=0)%s;',s, s));
                        newdb.(fds{j}) = permute(newdb.(fds{j}), backorder);
                        if ifsparse
                            newdb.(fds{j}) = sparse(newdb.(fds{j}));
                        end
                    end
                end
            end
        end            
        
        %append all fields as columns and return a table (cell or numeric
        %matrix)
        function table = fieldsToTable(datastruct, fields)            
            if any(~ismember(fields, fieldnames(datastruct)))
                error('some fields are not in datastruct');
            end
            nfd = length(fields);
            n = size( datastruct.(fields{1}), 1);
            ncol = zeros(1, n);
            isCellData = false(nfd, 1);
            for i = 1:nfd
                [n0, ncol(i)] = size( datastruct.(fields{i}));
                if n ~= n0
                    error('%s has different number of rows', fields{i});
                end                
                isCellData(i) = iscell( datastruct.(fields{i}) );                
            end
            
            if all(~isCellData)
                table = NaN(n, sum(ncol));
            else
                table = cell(n, sum(ncol));
            end
            isMixData = sum(isCellData) > 0 & sum(isCellData) < nfd;
            curCol = 0;
            for i = 1:nfd
                if isMixData && ~isCellData(i)
                    table(:, curCol+1:curCol+ncol(i)) = numarray2strarray(datastruct.(fields{i}));                    
                else
                    table(:, curCol+1:curCol+ncol(i)) = datastruct.(fields{i});
                end
                curCol = curCol + ncol(i);
            end
        end
        
        function newvar = allocateData(datatype, allocSize)
            if strcmp(datatype, 'cell')
                newvar = cell(allocSize);
                newvar(:) = {'NA'};
            else
                newvar = NaN(allocSize);
            end
        end
    end
    
end