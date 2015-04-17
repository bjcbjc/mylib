classdef GENOMEFUNC < handle
    
    %datastruct: NCBI
    methods (Static)        
        function [tf, idx] = overlapRegion(base, query)
            %base and query are matrices of three columns: chr, start, end
            %search if query is overlapped in base (tf=true or false);
            %idx are the indices of base that overlaps query
            %
            %query is a vector
            search = base(:,1) == query(1,1) & ...
                base(:,2)<=query(1,3) & query(1,2)<=base(:,3);
            tf = any(search);
            idx = find(search);
        end
        
        function tf = overlapRegionPerChrm(base, query)
            %base and query are matrices of two coloumn: start and end
            %
            %tf is a matrix of #base_row x #query_row
            %Assume base and query are on the same chromsome
            %
            if isempty(base) || isempty(query)
                tf = [];
                return
            end
            tf = bsxfun(@le, base(:,1), query(:,2)');
            tf = tf & bsxfun(@ge, base(:,2), query(:,1)');
        end
        
        function [tf, idx] = isWithinRegion(query, searchRegion)
            % query: [chr, start, (end)]
            % searchRegion: matrix, [chr, start, end]
            % tf: true or false column vector, same number of rows as query
            % idx:
            %   if query is a single locus, idx returns indices of all rows
            %   in searchRegion that contain query
            %   if query is a matrix, idx returns the same number of rows
            %   as query; each element contains the index of the first row
            %   in searchRegion that contains the row of the query
            %
            [nloc, locdim] = size(query);
            dim3 = min(3, locdim);
            if nloc == 1
                idx = query(:,1) == searchRegion(:,1) & ...
                    query(:,2) >= searchRegion(:,2) & ...
                    query(:,dim3) <= searchRegion(:,3);                            
                tf = any(idx);
                if nargout == 2                    
                    idx = find(idx);
                end
            else
                tf = false(nloc, 1);
                if nargout == 2
                    idx = zeros(nloc, 1);
                end
                for chrm = 1:25
                    subidx1 = query(:,1) == chrm;
                    subidx2 = searchRegion(:,1) == chrm;
                    mtx = bsxfun(@ge, query(subidx1,2), searchRegion(subidx2,2)') & ...
                        bsxfun(@le, query(subidx1,dim3), searchRegion(subidx2,3)');
                    tf(subidx1) = any(mtx, 2);
                    if nargout == 2
                        subidx1 = find(subidx1);
                        subidx2 = find(subidx2);
                        idx(subidx1(tf(subidx1))) = subidx2(arrayfun(@(x) find(mtx(x,:),1,'first'), find(tf(subidx1))));
                    end
                end                
            end
        end
        
        function [name, idx] = entrez2Name(datastruct, gid)
            newid = GENOMEFUNC.replaceOldID(datastruct, gid);
            [~, idx] = ismember(newid, datastruct.id);
            name = cell(length(gid), 1);
            name(idx~=0) = datastruct.name(idx(idx~=0));
            %for old id (retired), take the old name
            i = find(idx==0);
            [~, j] = ismember(newid(i), datastruct.history.oldid);
            name(i(j~=0)) = datastruct.history.name(j(j~=0));
            name(i(j==0)) = arrayfun(@num2str, newid(i(j==0)), 'uniformoutput', false);
        end
        
        function [gid, idx] = name2Entrez(datastruct, name, chrm, region)
            if nargin < 3, chrm = []; end            
            if nargin < 4, region = []; end
            if ischar(name)
                name = {name};
            end
            
            name = upper(name);            
            if isempty(chrm)
                fprintf('name mapping might not be unique\n');
                [~, idx] = ismember(name, upper(datastruct.name));
                notfound = find(idx == 0);
                if isfield(datastruct, 'alias')
                    [~, j] = ismember(name(idx==0), upper(datastruct.alias));
                    [~, gididx] = ismember(datastruct.alias_id(j(j~=0)), datastruct.id);
                    idx(notfound(j~=0)) = gididx;
                end
                notfound = find(idx == 0);
                if isfield(datastruct, 'locustag')
                    [~, j] = ismember(name(idx==0), upper(datastruct.locustag));
                    [~, gididx] = ismember(datastruct.locustag_id(j(j~=0)), datastruct.id);
                    idx(notfound(j~=0)) = gididx;
                end
                gid = NaN(length(name), 1);
                gid(idx ~= 0) = datastruct.id(idx( idx ~= 0 ));
            else
                %use chrm and region to narrow down mapping
                if ischar(chrm)
                    chrm = {chrm};
                end
                if iscellstr(chrm)                    
                    chrm = strrep(chrm, 'X|Y', '23');
                    chrm = strrep(chrm, 'X', '23');
                    chrm = strrep(chrm, 'Y', '24');
                    chrm = strrep(chrm, 'MT', '25');
                    chrm = strrep(chrm, 'Un', 'NaN');
                    chrm = strrep(chrm, 'chr', '');
                elseif isnumeric(chrm)
                    chrm = arrayfun(@num2str, chrm, 'uniformoutput', false);
                end
                
                gid = NaN(length(name), 1);
                
                query = strcat(name, '_', chrm);                
%                 key = strcat(datastruct.name, '_', arrayfun(@num2str, datastruct.chrm, ...
%                     'uniformoutput', false));
                key = strcat(upper(datastruct.name), '_', numarray2strarray(datastruct.chrm));
                [~, idx] = ismember(query, key);
                gid(idx~=0) = datastruct.id(idx(idx~=0));
                
                extralookup = {'alias', 'locustag'};
                for lookupidx = 1:length(extralookup)
                    notfound = find(idx == 0);
                    key = strcat(upper(datastruct.(extralookup{lookupidx})), ...
                        '_', numarray2strarray(datastruct.([extralookup{lookupidx} '_chrm'])));                    
                    [~, ii] = ismember(query(idx==0), key);
                    gid(notfound(ii~=0)) = datastruct.([extralookup{lookupidx}, '_id'])(ii(ii~=0));
                    [~, j] = ismember(datastruct.([extralookup{lookupidx}, '_id'])(ii(ii~=0)), datastruct.id);
                    idx(notfound(ii~=0)) = j;
                end                
                idx(idx == 0) = NaN;
            end
        end
        
        function newid = replaceOldID(datastruct, gid)
            i = find(~ismember(gid, datastruct.id));
            [~, j] = ismember(gid(i), datastruct.history.oldid);
            newid = gid;
            newid(i(j~=0)) = datastruct.history.newid(j(j~=0));
            %if a gene is retired, keep the old id
            newid(isnan(newid)) = gid(isnan(newid));
        end
        
        function [gid, idx] = closestGenes(datastruct, loc, whichend)
            %loc: [chrm, position1, position2 (optional)]
            %whichend: {'both' (default) , '5', '3'}
            %
            if nargin < 3, whichend = 'both'; end
            if size(loc,2) == 3
                loc(:,2) = nanmean(loc(:, 2:3), 2);
                loc(:,3) = [];
            end
            loc = double(loc);
            
            idx = NaN(size(loc,1), 1);
            uc = nanunique(loc(:,1));
            
            validloc = all(~isnan(loc), 2);
            
            if strcmpi(whichend, 'both')
                for ci = 1:length(uc)
                    chrmidx = find(datastruct.chrm == uc(ci));
                    dist1 = abs(bsxfun(@minus, datastruct.start(chrmidx), ...
                        loc(loc(:,1) == uc(ci) & validloc, 2)'));
                    dist2 = abs(bsxfun(@minus, datastruct.end(chrmidx), ...
                        loc(loc(:,1) == uc(ci) & validloc, 2)'));
                    dist1 = nanmin(dist1, dist2);
                    [~, mi] = nanmin(dist1);
                    idx(loc(:,1) == uc(ci) & validloc) = chrmidx(mi);
                end                
            elseif strcmpi(whichend, '5')
                for ci = 1:length(uc)
                    chrmidx = find(datastruct.chrm == uc(ci) & ~strcmp(datastruct.strand, ''));
                    dbloc = NaN(length(chrmidx), 1);
                    wstrand = strcmp(datastruct.strand(chrmidx), '+');
                    dbloc(wstrand) = datastruct.start(chrmidx(wstrand));
                    dbloc(~wstrand) = datastruct.end(chrmidx(~wstrand));
                    [~, mi] = nanmin( abs( bsxfun(@minus, dbloc, ...
                        loc(loc(:,1) == uc(ci) & validloc, 2)' ) ) );
                    idx(loc(:,1) == uc(ci) & validloc) = chrmidx(mi);
                end
            elseif strcmpi(whichend, '3')
                for ci = 1:length(uc)
                    chrmidx = find(datastruct.chrm == uc(ci) & ~strcmp(datastruct.strand, ''));
                    dbloc = NaN(length(chrmidx), 1);
                    wstrand = strcmp(datastruct.strand(chrmidx), '+');
                    dbloc(wstrand) = datastruct.end(chrmidx(wstrand));
                    dbloc(~wstrand) = datastruct.start(chrmidx(~wstrand));
                    [~, mi] = nanmin( abs( bsxfun(@minus, dbloc, ...
                        loc(loc(:,1) == uc(ci) & validloc, 2)' ) ) );
                    idx(loc(:,1) == uc(ci) & validloc) = chrmidx(mi);
                end
            else
                error('unknown specification of which end: %s', whichend);
            end
            
            gid = NaN(size(loc,2), 1);
            gid(~isnan(idx)) = datastruct.id(idx(~isnan(idx)));
        end
        
        function mapped = nearByGenes(datastruct, loc, whichend, window)
            %loc: [chrm, position1, position2 (optional)]
            %whichend: {'both' (default) , '5', '3'}
            %
            if nargin < 4, window = 2000; end
            if nargin < 3, whichend = 'both'; end
            if size(loc,2) == 3
                loc(:,2) = nanmean(loc(:, 2:3), 2);
                loc(:,3) = [];
            end
            loc = double(loc);
            
            mapped.gid = datastruct.id;
            mapped.nearby = sparse(length(datastruct.id), size(loc,1));
            uc = nanunique(loc(:,1));
            
            validloc = all(~isnan(loc), 2);
            
            if strcmpi(whichend, 'both')
                for ci = 1:length(uc)
                    chrmidx = find(datastruct.chrm == uc(ci));
                    locidx = find(loc(:,1) == uc(ci) & validloc)';
                    dist1 = abs(bsxfun(@minus, datastruct.start(chrmidx), ...
                        loc(locidx, 2)'));
                    dist2 = abs(bsxfun(@minus, datastruct.end(chrmidx), ...
                        loc(locidx, 2)'));
                    dist1 = nanmin(dist1, dist2);
                    for loci = 1:length(locidx)                        
                        mapped.nearby(chrmidx(dist1(:, loci) <= window), locidx(loci)) = 3;
                    end                    
                    [dist1, mi] = min(dist1);
                    coli = locidx(dist1 <= window); %closest and within window
                    mapped.nearby( sub2ind(size(mapped.nearby), mi(dist1 <= window), coli) ) = 1;
                    coli = locidx(dist1 > window); %closest but outside window
                    mapped.nearby( sub2ind(size(mapped.nearby), mi(dist1 > window), coli) ) = 2;
                end                
            elseif strcmpi(whichend, '5')
                for ci = 1:length(uc)
                    chrmidx = find(datastruct.chrm == uc(ci) & ~strcmp(datastruct.strand, ''));
                    locidx = find(loc(:,1) == uc(ci) & validloc)';
                    dbloc = NaN(length(chrmidx), 1);
                    wstrand = strcmp(datastruct.strand(chrmidx), '+');
                    dbloc(wstrand) = datastruct.start(chrmidx(wstrand));
                    dbloc(~wstrand) = datastruct.end(chrmidx(~wstrand));
                    dist = abs( bsxfun(@minus, dbloc, loc(locidx, 2)' ) );
                    for loci = 1:length(locidx)                        
                        mapped.nearby(chrmidx(dist(:, loci) <= window), locidx(loci)) = 3;
                    end                                        
                    [dist, mi] = min(dist);
                    coli = locidx(dist <= window); %closest and within window
                    mapped.nearby( sub2ind(size(mapped.nearby), mi(dist <= window), coli) ) = 1;
                    coli = locidx(dist > window); %closest but outside window
                    mapped.nearby( sub2ind(size(mapped.nearby), mi(dist > window), coli) ) = 2;
                end
            elseif strcmpi(whichend, '3')
                for ci = 1:length(uc)
                    chrmidx = find(datastruct.chrm == uc(ci) & ~strcmp(datastruct.strand, ''));
                    locidx = find(loc(:,1) == uc(ci) & validloc)';
                    dbloc = NaN(length(chrmidx), 1);
                    wstrand = strcmp(datastruct.strand(chrmidx), '+');
                    dbloc(wstrand) = datastruct.end(chrmidx(wstrand));
                    dbloc(~wstrand) = datastruct.start(chrmidx(~wstrand));
                    dist = abs( bsxfun(@minus, dbloc, loc(locidx, 2)' ) );
                    for loci = 1:length(locidx)                        
                        mapped.nearby(chrmidx(dist(:, loci) <= window), locidx(loci)) = 3;
                    end                 
                    [dist, mi] = min(dist);
                    coli = locidx(dist <= window); %closest and within window
                    mapped.nearby( sub2ind(size(mapped.nearby), mi(dist <= window), coli) ) = 1;
                    coli = locidx(dist > window); %closest but outside window
                    mapped.nearby( sub2ind(size(mapped.nearby), mi(dist > window), coli) ) = 2;
                end
            else
                error('unknown specification of which end: %s', whichend);
            end
            mapped.note = sprintf('1: closest and within %dbp;\n2: closest but outside %dbp;\n3: genes (except coded 1) that are within %dbp', ...
                window, window, window);
        end
        
        function q = idQuery(datastruct, query, returnFields, idname)
            if nargin < 4, idname = 'id'; end
            [~, i] = ismember(query, datastruct.(idname));
            n = length(query);
            if ~iscellstr(returnFields)            
                returnFields = {returnFields};
            end
            if any(cellfun(@(x) strcmp(x, 'all'), returnFields))
                returnFields = {'id', 'name', 'chrm', 'start', 'end', 'strand'};
            end
            
            q = cell(n, length(returnFields));
            allNumeric = true;
            for ki = 1:length(returnFields)
                if iscell(datastruct.(returnFields{ki}))
                    q(i==0, ki) = {'NotFound'};
                    if any(i ~= 0)
                        q(i~=0, ki) = datastruct.(returnFields{ki})(i(i~=0));
                    end
                    allNumeric = false;
                else
                    q(i==0, ki) = {NaN};
                    if any(i ~= 0)
                        q(i~=0, ki) = num2cell(datastruct.(returnFields{ki})(i(i~=0)));
                    end
                end                
            end  
            if allNumeric            
                q = cell2mat(q);
            end
        end        
        
        function q = nameQuery(datastruct, searchName, key)
            %name query might not work due to lack of uniqueness in gene
            %symbols            
            if length(unique(lower(datastruct.name))) < length(datastruct.name)
                fprintf('Gene symbols are not unique; some mapping may be incorrect\n');
            end
            if ~iscellstr(searchName)
                searchName = {searchName};
            end
            if ~iscellstr(key)
                key = {key};
            end
            n = length(searchName);
            [~,i] = ismember(lower(searchName), lower(datastruct.name));
            
            q = cell(n, length(key));
            allNumeric = true;
            for ki = 1:length(key)
                if iscell(datastruct.(key{ki}))
                    q(i==0, ki) = {'NotFound'};
                    q(i~=0, ki) = datastruct.(key{ki})(i(i~=0));
                    allNumeric = false;
                else
                    q(i==0, ki) = {NaN};
                    q(i~=0, ki) = num2cell(datastruct.(key{ki})(i(i~=0)));
                end
            end     
            if allNumeric
                q = cell2mat(q);
            end            
        end
        
        function g = genesInWindow(datastruct, chrm, sloc, eloc, key)
            id = datastruct.id(datastruct.chrm==chrm & ...
                datastruct.start <= eloc & ...
                datastruct.end >= sloc);
            if ischar(key)
                if strcmp(key, 'id')
                    g = id;
                end
            elseif iscell(key)
                if strcmp(key{1}, 'id') && length(key) == 1
                    g = id;
                end
            else
                g = GENOMEFUNC.idQuery(datastruct, id, key);
            end
        end
        
        function g = genesAround(datastruct, id, key, window)
            %for single gene 
            if nargin < 4,
                window = 10e6;
            end
            %loc = idquery(id, {'chrom', 'start', 'end'});
            i = datastruct.id == id;
            if ~any(i)
                g = [];
                return
            end
            loc = [datastruct.chrm(i) datastruct.start(i) datastruct.end(i)];
            g = GENOMEFUNC.genesInWindow(datastruct, loc(1), loc(2)-window, loc(3)+window, key);
        end
        
        function [tf, idx] = ismemberGencode(gencodeId1, gencodeId2, removedot)
            if nargin < 3, removedot = true; end
            if removedot
                gencodeId1 = regexp(gencodeId1, '\w+', 'match', 'once');
                gencodeId2 = regexp(gencodeId2, '\w+', 'match', 'once');
            end
            [tf, idx] = ismember(gencodeId1, gencodeId2);
        end
        
        function [overlap, idx1, idx2] = intersectGencode(gencodeId1, gencodeId2, removedot)
            if nargin < 3, removedot = true; end
            if removedot
                gencodeId1 = regexp(gencodeId1, '\w+', 'match', 'once');
                gencodeId2 = regexp(gencodeId2, '\w+', 'match', 'once');
            end
            [overlap, idx1, idx2] = intersect(gencodeId1, gencodeId2);
        end
        
        function result = searchAllOverlapGTF( loc, gtfFn, varargin )
            para.padChr = true;
            para.mt = 'M'; %or 'MT'
            para.bedtools = 'bedtools ';
            para.geneIdFormat = 'ENSG[\w\.]+';
            para.getUnique = false;
            para = assignpara(para, varargin{:});
            
            singleBase = false;
            if size(loc, 2) == 2
                loc(:,3) = loc(:,2);
                singleBase = true;
            end            
            randfn = sprintf('/tmp/rand%05d.bed',randi(10000,1));
            nTry = 1;
            while exist(randfn, 'file') ~= 0
                if nTry > 10
                    error('more than 10 tries for random file generation');                    
                end
                randfn = sprintf('/tmp/rand%05d.bed',randi(10000,1));
                nTry = nTry + 1;
            end
                
            loctxt = numarray2strarray(loc);
            loctxt(:,1) = strrep(loctxt(:,1), '23', 'X');
            loctxt(:,1) = strrep(loctxt(:,1), '24', 'Y');
            if strcmp(para.mt, 'M')
                loctxt(:,1) = strrep(loctxt(:,1), 'MT', 'M');
            elseif strcmp(para.mt, 'MT')
                loctxt(:,1) = strrep(loctxt(:,1), 'M', 'MT');
            end
            if para.padChr
                loctxt(:,1) = strcat('chr', loctxt(:,1));
            end
            tabwrite( randfn, loctxt);
            
            cmd = sprintf('%s intersect -wa -wb -a %s -b <(cut -f 1,4,5,9 %s)', ...
                para.bedtools, randfn, gtfFn);
            [status, output] = system(cmd);
            if status ~= 0                
%                 system(sprintf('rm -f %s',randfn));
                error('bedtools error %s, %d', randfn, status);                
            else
                t = textscan(output, '%s %s %s %*s %*s %*s %s', 'delimiter', '\t');
%                 lockey = strcat(loctxt(:,1), '-', loctxt(:,2), '-', loctxt(:,3));
                chrm = numericchrm( t{1} );
                result.locidx_start = gloc2index(chrm, str2double_fast(t{2}));                
                if ~singleBase
                    result.locidx_end = gloc2index(chrm, str2double_fast(t{3}));
                end
                result.geneId = regexp(t{4}, para.geneIdFormat, 'match', 'once');                
                system(sprintf('rm -f %s',randfn));
                %get unique entries
                if para.getUnique
                    if ~singleBase
                        [~, uidx] = unique(strcat(numarray2strarray(result.locidx_start), '_', ...
                            numarray2strarray(result.locidx_end), '_', ...
                            result.geneId));                        
                        result.locidx_end = result.locidx_end(udix);
                    else
                        [~, uidx] = unique(strcat(numarray2strarray(result.locidx_start), '_', ...                         
                            result.geneId));                        
                    end
                    result.locidx_start = result.locidx_start(uidx);                        
                    result.geneId = result.geneId(uidx);
                end
            end            
        end
        
    end
end