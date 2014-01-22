classdef AlleleCountData < handle
    
    properties
        filename;
        loc;
        numreads;
        ref;
        %indelcount;
        %indelseq;
        indels;
        acount;
        ntlabels = {'>', 'A', 'T', 'C', 'G', 'N', '*'; '<', 'a', 't', 'c', 'g', 'n', 'dummy'};
        
        patient;
        sample;
        minRead;
        maxDepth;
        minMapQ;
        stranded;
    end
    
    methods
        
        function obj = AlleleCountData(filename)
            obj.filename = filename;
            t = parseText(filename, 'ncolname', 0, 'nrowname', 0, 'numericcol', [2, 4]);
            obj.loc = [ str2double( ...
                strrep(strrep(strrep(strrep(strrep( t.text(:,1), 'chr', ''), 'X', '23'), 'Y', '24'), 'MT', '25'), 'M', '25')), ...
                t.numtext(:,1) ];
            obj.numreads = t.numtext(:,2);
            obj.ref = t.text(:,2);
            nloc = size(obj.loc, 1);
            %obj.indelseq = cell(nloc,1);
            %obj.indelcount = zeros(nloc, 2);
            obj.acount = zeros([nloc, size(obj.ntlabels')]);
            acountSize = size(obj.acount);
            query = regexp(t.text(:,3), '([\+-])(\d+)((??[>ATCGN<atcgn]{$2}))(\d+)', 'tokens');
            obj.indels = cellfun(@(x) vertcat(x{:}), query, 'unif', 0);
            query = regexprep(t.text(:,3), '([\+-])(\d+)((??[>ATCGN<atcgn]{$2}))(\d+)', '');
            query = regexp(query, '([>ATCGN<atcgn\*]+)(\d+)', 'tokens');
            query = cellfun(@(x) vertcat(x{:}), query, 'unif', 0);
%             for i = 1:nloc
%                 [~, idx] = ismember(query{i}(:,1), obj.ntlabels);
%                 [idx3, idx2] = ind2sub(size(obj.ntlabels), idx);
%                 obj.acount( sub2ind(acountSize, repmat(i, length(idx3), 1), idx2, idx3) ) = str2double( query{i}(:,2) );
%             end
            for idx2 = 1:size(obj.acount, 2)
                for idx3 = 1:size(obj.acount, 3)
                    validx = cellfun(@(x) find(strcmp(x(:,1), obj.ntlabels{idx3, idx2})), query, 'unif', 0);
                    idx1 = find(~cellfun(@isempty, validx));
                    if isempty(idx1), continue; end
                    val = str2double( arrayfun(@(a,b) query{a}{b,2}, idx1, cell2mat(validx(idx1)), 'unif', 0) );
                    n = length(idx1);
                    obj.acount( sub2ind(acountSize, idx1, ones(n,1)*idx2, ones(n,1)*idx3) ) = val;
                end
            end
        end
        function matfilename = addparainfo(obj)
            fn = regexp(obj.filename, '\S+/(\S)+$', 'tokens', 'once');
            fn = fn{1};
            if ~isempty(strfind(fn, 'Sample'))
                tokens = regexp(fn, 'Sample_([\w\-]+)_\S+.r(\d+).d(\w+).q(\w+).s(\w)', 'tokens', 'once');
            else
                tokens = regexp(fn, '([\w\-]+).\S+.r(\d+).d(\w+).q(\w+).s(\w)', 'tokens', 'once');
            end
            obj.sample = tokens{1};
            tmp = textscan(tokens{1}, '%s', 'delimiter', '-');
            obj.patient = tmp{1}{1};
            obj.minRead = str2double(tokens{2});
            obj.maxDepth = str2double(strrep(tokens{3}, 'k', '000'));
            obj.minMapQ = str2double(tokens{4});
            if strcmpi(tokens{5}, 't')
                obj.stranded = true;
            else
                obj.stranded = false;
            end
            matfilename = sprintf('acount.%s.r%s.d%s.q%s.mat', ...
                obj.sample, tokens{2}, tokens{3}, tokens{4} );
%             save(['data/' matfilename], 'aCountData');
%             clear aCountData
        end
        
    end
    
    methods (Static)
        function data = readTableFormatOutput(filename)
            f = fopen(filename, 'r');
            line = fgetl(f);
            fclose(f);
            token = textscan(line, '%s');
            ncol = length(token{1});
            t = parseText(filename, 'nrowname', 0, 'ncolname', 1, 'numericcol', [2, 4:ncol-3]);
            data.locidx = gloc2index( [ numericchrm(t.text(:,strcmp(t.colname, 'chr'))) t.numtext(:,strcmp(t.numcolname, 'pos')) ]);
            data.ref = t.text(:, strcmp(t.colname, 'ref'));
            data.numread = t.numtext(:, strcmp(t.numcolname, '#read'));
            ntbaseidx = ~ismember(t.numcolname, {'pos', '#read'});
            data.ntbase = t.numcolname( ntbaseidx);
            data.count = t.numtext(:, ntbaseidx);
            indelcol = find(strcmp(t.colname, 'indel'));
            data.indel = t.text(:, indelcol);
            data.indelcount = t.text(:, indelcol+1:indelcol+2);
            [data.indelcount, data.indel] = AlleleCountData.parseIndel(data);
        end
        
        function [count, indel] = indelCount(data, stranded, countmode)
            % stranded: true or false; def=false; count indel with strands
            %   or not; stranded only matters if the count table was
            %   generated with strand information; if false, indels are
            %   collapsed to non-stranded count before counting
            % countmode: {'all', 'max'}; 'all' sums up the counts of all
            %   indels; 'max' returns the maximum count of the indels;
            %   def='all'
            %
            % count: vector of indel counts
            % indel: only returns values if 'max' is used as countmode,
            %   indicating which indel has the max count
            %
            if nargin < 3, countmode = 'all'; end
            if nargin < 2, stranded = false; end
            if strcmpi(countmode, 'all')
                indel = [];
                count = cellfun(@sum, data.indelcount);                
                if ~stranded 
                    count = sum(count,2);
                end
            elseif strcmpi(countmode, 'max')
                indel = cell(length(data.ref), 1);
                count = zeros(length(data.ref), 1);
                indel(:) = {''};
                if ~stranded && size(data.indelcount,2) > 1
                    tmp = cellfun(@(x,y) x+y, data.indelcount(:,1), data.indelcount(:,2), 'unif', 0);
                else
                    tmp = data.indelcount;
                end
                idx = cellfun(@length, data.indel) > 0;
                [count(idx), midx] = cellfun(@max, tmp);
                idx = find(idx);
                indel(idx) = arrayfun(@(i) data.indel{idx(i)}{midx(i)}, 1:length(idx), 'unif', 0);                
            else
                error('Unknown count mode %s', countmode);
            end
        end
        
        function [count, maxAllele] = maxNonRef(data, stranded)
            if nargin < 2, stranded = false; end
            if ~stranded
                data = AlleleCountData.collapseStrand(data);
            end            
            colidx = ismember(data.ntbase, {'A','T','C','G','a','t','c','g','*'});
            ntbase =data.ntbase(colidx);
            tmp = data.count(:, colidx);
            [~, refidx] = ismember(data.ref, ntbase);
            tmp( sub2ind(size(tmp), (1:length(refidx))', refidx) ) = 0;
            [count, altidx] = max(tmp, [], 2);
            maxAllele = ntbase( altidx);
        end
        
        function dataNoStrand = collapseStrand(data)
            dataNoStrand = data;
            if all(ismember({'A','a'}, data.ntbase))
                [~, keep] = ismember({'>','A','T','C','G','N'}, data.ntbase);
                [~, remove] = ismember({'<', 'a','t','c','g','n'}, data.ntbase);
                dataNoStrand.count(:, keep) = dataNoStrand.count(:, keep) + dataNoStrand.count(:, remove);
                dataNoStrand.ntbase(remove) = [];
                dataNoStrand.count(:, remove) = [];
                dataNoStrand.ref = upper(dataNoStrand.ref);
            end
            if size(data.indelcount, 2) > 1
                dataNoStrand.indelcount = cellfun(@(x,y) x+y, data.indelcount(:,1), data.indelcount(:,2), 'unif', 0);
            end            
        end
        
        function [count, indel] = parseIndel(data)
            n = length(data.ref);
            count = cell(n, size(data.indelcount,2));
            indel = cell(n, 1);
            idx = ~strcmp(data.indel, '');
            indel(~idx) = {''};
            count(~idx, :) = {0};
            indel(idx) = regexp(data.indel(idx), ',', 'split');
            count(idx, :) = cellfun(@(x) str2double_fast(x), ...
                regexp( data.indelcount(idx, :), ',', 'split'), 'unif', 0);
        end
    end
end

