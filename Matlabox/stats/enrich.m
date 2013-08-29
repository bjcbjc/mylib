function [pstat cat] = enrich(data, orf, method, nperm, geneset, pcut)
    %data can be binary matrix or real-value matrix
    %
    %Input:
    %   data: #genes * #series (#conditions); cell array of orfs
    %   orf: all ORFs in data, a vector; length should be the same as #genes
    %       in data if data is a matrix; if omitted, assume to the data has the
    %       same number of ORFs as that in geneset
    %
    %   method: {'hype': hypergeometric, 'KS': Kolmogorov-Smirnov, 'GSEA':
    %       GSEA-like (since we don't have phenotypes)
    %   nperm: #permutation for GSEA-like only
    %   geneset: cell array of strings {'GO','KEGG',...}
    %   pcut: cutoff for pval, def = 0.001;
    %
    %return: 
    %   pstat: matrix, #passed test x 4 [condi pval X Y], X is the number of
    %       genes in the condition and in the category, and Y is the number 
    %       of genes in the category (and in the data)
    %   cat: cell    
    %   
    %
    %revised: 06/24/08, take the number of overlap as the size of base
    %revised: 03/05/08, take all GO's genes
    %

    if nargin < 2, orf = []; end
    if nargin < 3, method = 'h'; end
    if nargin < 4, nperm = 1000; end
    if nargin < 5, geneset = {'go'}; end
    if nargin < 6, pcut = 0.001; end
      
    if ~iscell(geneset), geneset = {geneset}; end    
    
    datacopy = [];
    if ~iscell(data) && ~isempty(orf)
        if size(data,1) ~= length(orf)
            error('#orfs ~= size(data,1)\n');
        end
    elseif ~isempty(orf)
        [tmp i] = intersect(orf, data);
        data = zeros(length(orf),1);
        data(i) = 1;
    elseif iscell(data)
        datacopy = data;
    end
    
    %returned results
    pstat = zeros(0,4);
    cat = cell(0,1);
    
    for gseti = 1:length(geneset)
        switch lower(geneset{gseti})
            case 'go'
                s = load('C:\Documents and Settings\BJC.COLUMBIA-391490\Desktop\Yeast_growth\DATA\GO\GO.mat');
                baseset = s.GO;
            otherwise
                warning('unknown gene set %s\n',geneset{gseti});
                continue;
        end
        
        if isempty(orf)
            orf = baseset.orf;
            if iscell(datacopy)
                [tmp i] = intersect(orf, datacopy);
                data = zeros(length(orf),1);
                data(i) = 1;
            end
        end                
        
        [tmp basegi datagi] = intersect(baseset.orf, orf);
        
        %only take those categories containing more than 15 genes
        %cati = gset(GO.mtx(gi,:), 15);
        
        filtered_data = data(datagi,:);        
        [ngene, ncond] = size(filtered_data);
        if ngene == 1, error('#genes should > 1\n'); end
                
        if ~isempty(strfind(lower(method),'h'))
            [pval esize csize] = hygeEnrich(filtered_data, baseset.mtx(basegi,:));

        end
        
        [condi cati] = find(pval <= pcut);
        %cat = [cat; baseset.cat(cati)];
        for i = 1:length(condi)
            if csize(cati(i)) < 4
                continue
            end
            pstat(end+1,:) = [condi(i), pval(condi(i),cati(i)), ...
                esize(condi(i)), csize(cati(i))];
            cat{end+1,1} = baseset.cat{cati(i)};
        end
    end
    [pstat si] = sortrows(pstat,[2 1]);
    cat = cat(si);
    
end
%         elseif ~isempty(strfind(lower(method),'k'))        
%             pval = KStestMtx(dataInGO, GO.mtx(:,cati));
%         elseif ~isempty(strfind(lower(method),'g'))        
%             [pval, res] = gsealike(dataInGO, GO.mtx(:,cati),nperm);

function [p enrich_size cat_size] = hygeEnrich(data, cats)
    [ngene, ncond] = size(data);
    ncat = size(cats,2);
    p = NaN(ncond, ncat);
    cat_size = sum(cats);
    enrich_size = sum(data);
    for j = 1:ncat
        overlap = sum(data & repmat(cats(:,j),1,ncond));
        p(:,j) = 1 - hygecdf(full(overlap), ngene, cat_size(j), full(enrich_size));
    end
end

function [p res] = gsealike(data, cats, nperm)
    parap = 1;        
    [p res] = permTest(nperm, 155, 'tail1', @gseaMorph, data, cats, 1);
end

function [coli] = gset(mtx, minsize)
    %filter gset, return indices for selected categories (column indices)
    %minsize: min size of a gene set    
    [ngene ncat] = size(mtx);
    coli = [];
    for i = 1:ncat
        if length(find(mtx(:,i)==1)) >= minsize
            coli = [coli i];
        end
    end
end