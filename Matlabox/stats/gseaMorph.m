function [ES] = gseaMorph(data, cats, parap)

    if nargin < 3
        parap = 1;
    end    
    
    [ngene, ncond] = size(data);
    ncat = size(cats,2);
    ES = NaN(ncond, ncat);
        
    Nh = sum(cats); %size of genesets    
    [sdata, si] = sort(data,'descend');
    sdata = abs(sdata) .^ parap;
    for i = 1:ncond        
        tmpdata = repmat(sdata(:,i),1,ncat);
        Nr = sum(  tmpdata .* cats(si(:,i),:)  );
        hits = cumsum(  tmpdata .* cats(si(:,i),:)  );
        miss = cumsum( ~cats(si(:,i),:) ) ./ repmat((ngene - Nh),ngene,1);
        ESmtx = hits - miss;
        [ES(i,:) ESi] = max(abs(ESmtx));
        ES(i,:) = ESmtx(ESi + (0:ngene:(ncat-1)*ngene));
    end
end