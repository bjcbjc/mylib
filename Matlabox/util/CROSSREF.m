classdef CROSSREF < handle
    
    methods (Static = true)
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
                    result = CROSSREF.allocateData(class(data.(queryFields{1})), [n, dim2]);
                    result(idx~=0, :) = data.(queryFields{1})(idx(idx~=0 ), :);
                elseif keyDim == dim2
                    result = CROSSREF.allocateData(class(data.(queryFields{1})), [dim1, n]);
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
                    data.(newfd) = CROSSREF.allocateData(class(copyFrom.(copyFields{i})), [n, dim2]);
                    data.(newfd)(idx ~= 0, :) = copyFrom.(copyFields{i})( idx(idx~=0), : );
                elseif length(copyFrom.(copyKey)) == dim2
                    data.(newfd) = CROSSREF.allocateData(class(copyFrom.(copyFields{i})), [dim1, n]);
                    data.(newfd)(:, idx ~= 0) = copyFrom.(copyFields{i})( :, idx(idx~=0));
                else
                    error('cannot determine which dimension to match for %s', copyFields{i});
                end
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