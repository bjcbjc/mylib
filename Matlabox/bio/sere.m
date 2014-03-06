function score = sere2(readcount, minReadCount, pairwise)
    % simple error ratio estimate (Schulze et al BMC Genomics 2012)
    %
    % readcount: #gene x #sample
    % minReadCount: threshold to filter genes
    % pairwise: calculate score for each pair of samples, return #sample x
    %   #sample matrix; if false, then calculate the score across ALL
    %   samples
    
    if nargin < 2, minReadCount = 1; end
    if nargin < 3, pairwise = true; end
    
    if ~pairwise
        %total is calculated before filtering (why?)
        total = sum(sum(readcount));
        nReadAcrossSample = sum(readcount, 2);
        nReadAcrossGene = sum(readcount, 1);
        
        expectation = bsxfun(@times, nReadAcrossSample, nReadAcrossGene) ./ total;
        
        valid = nReadAcrossSample > minReadCount;
        nsample = size(readcount, 2);        
        %single estimate
        score = sqrt( sum(sum( (readcount(valid,:)-expectation(valid,:)).^2 ./ expectation(valid,:) )) ./ sum(valid) ./ (nsample-1));                
    elseif minReadCount < 0 %no filtering: much faster for large dataset
        [ngene, nsample] = size(readcount);
        nReadAcrossGene = sum(readcount, 1);
        score = zeros(nsample,nsample);
        sqcount = readcount.^2;
        for i = 1:nsample-1
            for j = i+1:nsample
                nReadAcrossSample = sum(readcount(:, [i, j]), 2);
                total = sum(nReadAcrossGene([i,j])); 
                crsprod = bsxfun(@times, nReadAcrossSample, nReadAcrossGene([i,j]));
                score(i,j) = sqrt(total * (sum(sum(sqcount(:,[i,j]) ./ crsprod)) -1) ./ ngene);
                %score(i,j) = sqrt( sum(sum( (readcount(valid,[i,j])-expectation(valid,:)).^2 ./ expectation(valid,:) )) ./ sum(valid) );
                score(j,i) = score(i,j);
            end
        end
%         for i = 1:nsample-1
%             nReadAcrossSample = bsxfun(@plus, readcount(:,i), readcount(:,i+1:end)); %ngene x test            
%             total = sum(nReadAcrossSample, 1);
%             ns = zeros(1, nsample-i, 2);
%             ns(:,:,1) = nReadAcrossGene(i);
%             ns(:,:,2) = nReadAcrossGene(i+1:end);            
%             crsprod = bsxfun(@times, nReadAcrossSample, ns);  %gene x test x 2
%             sqobs = zeros(ngene, nsample-i, 2);
%             sqobs(:,:,1) = repmat(sqcount(:,i),1, nsample-i);
%             sqobs(:,:,2) = sqcount(:,i+1:end);
%             score(i, i+1:end) = sqrt(total .* (sum(sum(sqobs./crsprod,3),1) -1) ./ ngene);
%             score(i+1:end, i) = score(i, i+1:end);
%         end
    else                
        [ngene, nsample] = size(readcount);
        nReadAcrossGene = sum(readcount, 1);
        score = zeros(nsample,nsample);
        sqcount = readcount.^2;
        for i = 1:nsample-1
            for j = i+1:nsample
                nReadAcrossSample = sum(readcount(:, [i, j]), 2);
                total = sum(nReadAcrossGene([i,j]));
                %expectation = bsxfun(@times, nReadAcrossSample, nReadAcrossGene([i,j])) ./ total;                
                valid = nReadAcrossSample > minReadCount;                                
                crsprod = bsxfun(@times, nReadAcrossSample(valid), nReadAcrossGene([i,j]));
                score(i,j) = sqrt(total * (sum(sum(sqcount(valid,[i,j]) ./ crsprod)) -1) ./ sum(valid));
                %score(i,j) = sqrt( sum(sum( (readcount(valid,[i,j])-expectation(valid,:)).^2 ./ expectation(valid,:) )) ./ sum(valid) );
                score(j,i) = score(i,j);
            end
        end
%         for i = 1:nsample-1
%             nReadAcrossSample = bsxfun(@plus, readcount(:,i), readcount(:,i+1:end)); %ngene x test
%             valid = nReadAcrossSample > minReadCount;
%             nReadAcrossSample( ~valid ) = 0;
%             ns = zeros(1, nsample-i, 2);
%             ns(:,:,1) = nReadAcrossGene(i);
%             ns(:,:,2) = nReadAcrossGene(i+1:end);
%             ns = bsxfun(@rdivide, ns, sum(ns, 2));            
%             expectation = bsxfun(@times, nReadAcrossSample, ns);  %gene x test x 2
%             obs = zeros(ngene, nsample-i, 2);
%             obs(:,:,1) = repmat(readcount(:,i),1, nsample-i);
%             obs(:,:,2) = readcount(:,i+1:end);
%             s = sum( (obs - expectation).^2 ./ expectation, 3); %gene x test
%             s(~valid) = 0;
%             score(i, i+1:end) = sqrt( sum(s,1) ./ sum(valid));
%             %score(i, i+1:end) = sqrt(shiftdim(sum(sum( (obs - expectation).^2 ./ expectation, 1), 2),1) ./sum(valid));
%             score(i+1:end, i) = score(i, i+1:end);
%         end        
    end
