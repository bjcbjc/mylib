classdef genequery
    
    
    properties
        geneloc
        datatype
    end
    
    methods
        function obj = genequery(genelocfn, datatype)
            if nargin < 2
                datatype = 'gene';
            end
            obj.datatype = datatype;
            if nargin < 1
                if strcmpi(datatype, 'gene')
                    tmp = load('/Users/bjc/Desktop/Projects/DATA/cancerdb/geneloc.mat');
                elseif strcmpi(datatype, 'mir')
                    tmp = load('/Users/bjc/Desktop/Projects/DATA/cancerdb/mirloc.mat');
                else
                    error('unknown data type %s', datatype);
                end
                fd = fieldnames(tmp);
                obj.geneloc = tmp.(fd{1});
            else
                if ischar(genelocfn)
                    tmp = load(genelocfn);
                    fd = fieldnames(tmp);
                    obj.geneloc = tmp.(fd{1});
                else
                    obj.geneloc = genelocfn;
                end
            end
        end
        
        function q = idquery(obj, id, key)
            [~, i] = ismember(id, obj.geneloc.id);
            n = length(id);
            if ~iscellstr(key)            
                key = {key};
            end
            if any(cellfun(@(x) strcmp(x, 'all'), key))
                if strcmpi(obj.datatype, 'gene')
                    key = {'id', 'name', 'loc', 'arm', 'chrom', 'start', 'end'};
                elseif strcmpi(obj.datatype, 'mir')
                    key = {'id', 'name', 'chrom', 'start', 'end', 'strand'};
                end
            end
            
            q = cell(n, length(key));
            for ki = 1:length(key)
                switch key{ki}                    
                    case {'name', 'loc', 'arm', 'strand'}
                        q(i==0, ki) = {'NotFound'};
                        q(i~=0, ki) = obj.geneloc.(key{ki})(i(i~=0));
                    case {'chrom', 'start', 'end', 'id'}
                        q(i==0, ki) = {NaN};
                        q(i~=0, ki) = num2cell(obj.geneloc.(key{ki})(i(i~=0)));                    
                end
            end  
            if strcmpi(obj.datatype, 'gene')
                if ~any(ismember({'name', 'loc', 'arm'}, key))
                    q = cell2mat(q);
                end
            elseif strcmpi(obj.datatype, 'mir')
                if ~any(ismember({'name','id', 'strand'}, key))
                    q = cell2mat(q);
                end
            end
        end        
        
        function q = namequery(obj, name, key)
            %name query might not work due to lack of uniqueness in gene
            %symbols
            if ~strcmpi(obj.datatype, 'gene')
                fprintf('only support gene name query');
                q = {};
                return
            end
            if length(unique(lower(obj.geneloc.name))) < length(obj.geneloc.name)
                fprintf('Gene symbols are not unique; some mapping may be incorrect\n');
            end
            if ~iscellstr(name)
                name = {name};
            end
            if ~iscellstr(key)
                key = {key};
            end
            n = length(name);
            [~,i] = ismember(lower(name), lower(obj.geneloc.name));
            
            q = cell(n, length(key));
            for ki = 1:length(key)
                switch key{ki}                    
                    case {'name', 'loc', 'arm'}
                        q(i==0, ki) = {'NotFound'};
                        q(i~=0, ki) = obj.geneloc.(key{ki})(i(i~=0));
                    case {'chrom', 'start', 'end', 'id'}
                        q(i==0, ki) = {NaN};
                        q(i~=0, ki) = num2cell(obj.geneloc.(key{ki})(i(i~=0)));                    
                end
            end     
            if ~any(ismember({'name', 'loc', 'arm'}, key))
                q = cell2mat(q);
            end            
        end
        
        function g = genesinwindow(obj, chrm, sloc, eloc, key)
            id = obj.geneloc.id(obj.geneloc.chrom==chrm & ...
                obj.geneloc.start <= eloc & ...
                obj.geneloc.end >= sloc);
            if ischar(key)
                if strcmp(key, 'id')
                    g = id;
                end
            elseif iscell(key)
                if strcmp(key{1}, 'id') && length(key) == 1
                    g = id;
                end
            else
                g = obj.idquery(id, key);
            end
        end
        
        function g = genesaround(obj, id, key, window)
            %for single gene 
            if nargin < 4,
                window = 10e6;
            end
            %loc = obj.idquery(id, {'chrom', 'start', 'end'});
            i = obj.geneloc.id == id;
            if ~any(i)
                g = [];
                return
            end
            loc = [obj.geneloc.chrom(i) obj.geneloc.start(i) obj.geneloc.end(i)];
            g = obj.genesinwindow(loc(1), loc(2)-window, loc(3)+window, key);
        end
        
    end
    
end

