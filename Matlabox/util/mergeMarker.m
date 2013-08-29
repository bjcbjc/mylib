function indexes = mergeMarker(genotype, thres, markerindexes, indexthres, markercell)
    %genotype: #sample x #marker
    %thres: threshold for correlation
    %markerindexes: indexes of these markers (optional)
    %markercell: markers (optional)
    %indexthres: thres of maximum diff of indexes of markers (optional)    
    %
    %indexes: #marker x 1, merged indexes, eg, if 1~3 being merged, they will have the
    %same index (eg 1)
    
    if nargin < 3, markerindexes = []; end
    if nargin < 4, indexthres = 10; end
    if nargin < 5, markercell = {}; end    
    
    if ~isempty(markercell)
        chrm = getchrm(markercell);
    else
        chrm = [];
    end
    
    %only look at adjacent markers; c has nmarker-1 elements
    %c = diag(corr(genotype),1); 
    c = corr(genotype);
    nmarker = size(genotype,2);
    indexes = ones(nmarker,1);
    mask = c>=thres;    
    cur = 1;
    for i = 2:nmarker
        %if mask(i-1) %correlated
        smi = find(indexes==cur,1,'first');
        if nnz(mask(smi:i, smi:i)) == numel(mask(smi:i, smi:i))
            if ~isempty(chrm) %check chrm
                if chrm(i) ~= chrm(i-1)
                    cur = cur + 1;
                    indexes(i) = cur;
                    continue;
                end
            end
            if ~isempty(markerindexes) %check indexes
                %if markerindexes(i) - markerindexes(i-1) > indexthres
                if markerindexes(i) - markerindexes(smi) > indexthres
                    cur = cur + 1;
                    indexes(i) = cur;
                    continue;
                end
            end                
            indexes(i) = cur;
        else
            cur = cur + 1;
            indexes(i) = cur;
        end
    end
end

function chrm = getchrm(markercell)
    nm = length(markercell);
    chrm = zeros(nm,1);
    for i = 1:nm
        chrm(i) = str2num(strtok(strrep(markercell{i},'M',''),'_'));
    end
end