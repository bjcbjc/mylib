classdef DATAFUNC < handle
    
    methods (Static)
        
        function outputstr = translateByTable(inputstr, oldstrlist, newstrlist)
            % inputstr: string or cell array of strings
            % oldstrlist: list of old string; n x 1 or, n x 2 cell arrays, with
            %   the first column as old strings and the second column as new
            %   strings
            % newstrlist: list of new string, must be of the same length as
            %   oldstrlist
            %            
            if nargin < 3, newstrlist = {}; end
            if ischar(oldstrlist) && ischar(newstrlist)
                outputstr = strrep(inputstr, oldstrlist, newstrlist);
                return
            end
            if isempty(newstrlist) && ~isvector(oldstrlist)
                newstrlist = oldstrlist(:,2);
                oldstrlist(:,2) = [];
            end
            if ~isvector(oldstrlist) || ~isvector(newstrlist)
                error('strings for lookups are not list');                
            end
            n = length(oldstrlist);
            if length(newstrlist) ~= n
                error('string lists for lookups are not of the same size');
            end
            outputstr = inputstr;
            for i = 1:n
                outputstr = strrep(outputstr, oldstrlist{i}, newstrlist{i});
            end
        end

        %take unique values of two columns and return binary matrix
        function res = data2mtx(data1, data2, varargin)
            if nargin < 2, data2 = []; end
            para.default = 0;
            para.sparse = false;
            para = assignpara(para, varargin{:});
            [~, ncol1] = size(data1);
            [~, ncol2] = size(data2);
            if ncol1 + ncol2 ~= 2
                error('input must have two columns of data');
            end
            if ncol1 == 2
                data2 = data1;
                data1 = data1(:,1);
            end
            [res.rowLabel, ~, rowIdx] = unique(data1);
            [res.colLabel, ~, colIdx] = unique(data2);
            nRow = length(res.rowLabel);
            nCol = length(res.colLabel);
            if para.sparse
                res.mtx = sparse(rowIdx, colIdx, 1, nRow, nCol);
            else                
                res.mtx = NaN(nRow, nCol);
                res.mtx(:) = para.default;
                res.mtx( sub2ind([nRow, nCol], rowIdx, colIdx) ) = 1;
            end
        end
        
        function [binmtx, label] = cat2bin(categoricaldata)
            if ~isvector(categoricaldata)
                error('input categorical data must be a vector');
            end
            
            n = length(categoricaldata);
            [label, ~, uidx] = unique(categoricaldata);
            nl = length(label);
            binmtx = false(n, nl);
            binmtx( sub2ind( [n, nl], (1:n)', uidx ) ) = true;
        end
        
        %return data(index), but fill in NaN or NA for index=0 or NaN index
        function res = getDataByIndex(data, index, indexDim)
            if nargin < 3, indexDim = 1; end            
            if indexDim ~= 1
                newOrder = [indexDim, setdiff(1:ndims(data), indexDim)];
                [~, revertOrder] = sort(newOrder);
                data = permute(data, newOrder);
            end
            valid = index ~= 0 & ~isnan(index);
            if all(valid)
                res = data(index, :);
            else
                n = length(index);
                resSize = size(data);
                resSize(1) = n;
                res = DATAFUNC.allocateData(class(data), resSize);
                res(valid,:) = data(index(valid),:);
            end
            if indexDim ~= 1
                %revert back
                res = permute(res, revertOrder);
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