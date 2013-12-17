classdef TCGAParser < handle
   
    methods (Static)
        function data = readsdrf(fn)
            data = parseText(fn, 'ncolname',1, 'nrowname', 0);
        end
        function fns = matchID(fns, sdrf)
            fns(:,2) = {''};
            if ischar(sdrf)
                sdrf = TCGAParser.readsdrf(sdrf);
            end
            tmp = sdrf.text;
            idcol = strcmp(sdrf.colname, 'Comment [TCGA Barcode]');
            if sum(idcol) ~= 1
                error('number of ID columns: %d', sum(idcol));
            end
            rmi = cellfun(@isempty, strfind(tmp(:,idcol), 'TCGA'));
            tmp(rmi, :) = [];
            [~, i] = ismember(fns(:,1), tmp(:));
            [rowidx, ~] = ind2sub(size(tmp), i(i~=0));
            fns(i~=0,2) = tmp(rowidx, idcol);
            if any(i==0)
                fprintf('miss some %d samples\n', sum(i==0));
            end
        end
        
        function data = reaExpArray(fns, fndir)
            if nargin < 2, fndir = ''; end
            fns(:,1) = strcat(fndir, fns(:,1));
            n = size(fns, 1);            
            for i = 1:n                
                t = parseText(fns{i,1}, 'nrowname', 1, 'ncolname', 2, 'numeric', true, 'treatasempty', {'na','NA','null'});
%                 if ~strcmp(t.colname{1}, fns{i,2})
%                     error('inconsist ID %s, %s', t.colname{1}, fns{i,2});
%                 end
                if i == 1
                    data.gid = t.rowname;
                    data.exp = NaN(length(data.gid), n);
                    data.exp(:,i) = t.text;
                else
                    if any(~strcmp(data.gid, t.rowname))
                        error('file %d has different gene names', i);
                    end
                    data.exp(:,i) = t.text;
                end
            end            
            data.sample = fns(:,2);
        end
        
        function data = readRnaSeq(fns, fndir)            
            if nargin < 2, fndir = ''; end
            fns(:,1) = strcat(fndir, fns(:,1));
            n = size(fns,1);
            data.gid = [];
            tmpdata = [];
            normcount = false(n,1);
            for i = 1:n
                if ~isempty(strfind(fns{i,1}, 'normalized'))
                    sampledata = TCGAParser.readRnaSeq_genes_normalized(fns{i,1});
                    normcount(i) = true;
                else
                    sampledata = TCGAParser.readRnaSeq_genes(fns{i,1});                    
                end
                if i == 1
                    data.gid = sampledata.gid;
                    tmpdata = NaN(length(data.gid), n);
                end
                if any(data.gid ~= sampledata.gid)
                    error('file %d has different gene ids', i);
                end
                tmpdata(:,i) = sampledata.exp;
            end
            sample1 = fns(normcount, 2);
            sample2 = fns(~normcount, 2);
            if any(~strcmp(sample1, sample2))
                error('not all samples have raw and normalized rna counts');
            end
            data.sample = sample1;
            data.rnaCount = tmpdata(:, ~normcount);
            data.rnaNormCount = tmpdata(:, normcount);
        end
        
        function data = readRnaSeq_genes_normalized(fn)  
            t = parseText(fn, 'ncolname', 1, 'nrowname', 1, 'numeric', true, 'treatasempty', {'na','NA','null'});
            data.gid = str2double_fast(cellfun(@(b) b{1}, regexp(t.rowname, '[\w\?]\|(\d+)', 'tokens', 'once'),'unif',0)); 
            data.exp = t.text;
        end
        function data = readRnaSeq_genes(fn)
            t = parseText(fn, 'ncolname', 1, 'nrowname', 0, 'numericcol',[2,3], 'treatasempty', {'na','NA','null'});
            data.gid = str2double_fast(cellfun(@(b) b{1}, regexp(t.text(:,1), '[\w\?]\|(\d+)', 'tokens', 'once'),'unif',0)); 
            data.exp = t.numtext(:, strcmp(t.numcolname, 'raw_count'));
        end
    end
end