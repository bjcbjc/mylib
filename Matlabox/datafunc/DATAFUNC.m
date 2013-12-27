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
    end
    
end