classdef MembershipDataset < dataset
        
    
    methods
        function obj = MembershipDataset(colname, membershipMtx, varargin)
            %varargin: {'names'}, matrix, {'names'}, cell_arrays          
            colname = regexprep(colname, '(\W)+', '_');
            obj@dataset({sparse(logical(membershipMtx)), colname{:}});
            
            UserData.anntNames = colname;
            if ~isempty(varargin) 
                for i = 1:2:length(varargin)
                    if iscell(varargin{i+1})
                        obj = horzcat(obj, cell2dataset(varargin{i+1}, 'varnames', varargin{i}) );
                    else
                        obj = horzcat(obj, mat2dataset(varargin{i+1}, 'varnames', varargin{i}) );
                    end
                end
            end
            UserData.entityNames = setdiff(obj.Properties.VarNames, UserData.anntNames);
            obj = set(obj, 'UserData', UserData);
        end
        
        function res = annt(obj, rowidx, names)
            %return cell arrays containing names that rowidx belongs to
            %if arg names is given, limit returned names to be in this
            %scope
            
            if nargin < 3
                names = obj.Properties.UserData.anntNames;
            end
            s.type = '()';
            s.subs = {rowidx, names};
            mem = logical(double(subsref(obj, s)));
            res = arrayfun(@(i) names( mem(i,:) ), (1:size(mem,1))', 'unif', 0);            
        end
        
        function [varargout] = subsref(obj, S)
            if strcmp(S(1).type, '.')
                if ismethod(obj, S(1).subs)
                    [varargout{1:nargout}] = eval(sprintf('%s(obj,S(2).subs{:});', S(1).subs) );
                    return
                elseif ismember(S(1).subs, obj.Properties.UserData.entityNames) && length(S) > 1
                    idx = any(logical(double(subsref@dataset(obj, S(2:end)))), 2);
                    S(2:end) = [];
                    S(2).type = '()';
                    S(2).subs = {idx};
                    [varargout{1:nargout}] = subsref@dataset(obj, S);            
                    %varargout{1} = obj.(S(1).subs)();
                    return
                end
            end
            [varargout{1:nargout}] = subsref@dataset(obj, S);            
        end
        
    end
end