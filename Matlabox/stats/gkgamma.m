function gamma = gkgamma(ordinalData, ordinalData2)
    %ordinalData: #obs x #variable
    %   ordinalData are numeric discrete data but the numeric values
    %   suggest quantitative relationships
    %

    
    
    %create count tables; don't use crosstab; it's slow
    ulabel = unique(ordinalData(:));
    [sampleSize, nvar] = size(ordinalData);
    [sampleSize2, nvar2] = size(ordinalData2);
    if sampleSize ~= sampleSize2
        error('sample size not the same');
    end
    n = length(ulabel);
    table = zeros(n*n, nvar, nvar2);
    for vi = 1:nvar
        for i = 1:length(ulabel)
            for j = 1:length(ulabel)
                %count(i,j) for vi variable and others
                tableRowIdx = sub2ind([n,n],i,j);
                table(tableRowIdx, vi,:) = sum(bsxfun(@and, ordinalData(:,vi)==ulabel(i), ordinalData2==ulabel(j)), 1);
                %table(tableRowIdx, vi,vi:end) = sum(bsxfun(@and, ordinalData(:,vi)==ulabel(i), ordinalData(:,vi:end)==ulabel(j)), 1);
                %table(tableRowIdx, vi:end, vi) = table(tableRowIdx, vi, vi:end);
            end
        end
    end
    
    %calculate concordance and disconcordance
    con = zeros(size(table));
    dis = zeros(size(table));
    dummyidx = reshape(1:n*n, n, n);
    for tableRowIdx = 1:n*n
        [ui, uj] = ind2sub([n,n], tableRowIdx);
        idx = [reshape(dummyidx( 1:ui-1, 1:uj-1), (ui-1)*(uj-1),1); ...
            reshape(dummyidx(ui+1:end, uj+1:end), (n-ui)*(n-uj), 1)];
        con(tableRowIdx, :, :) = sum( table(idx, :,:), 1);
        idx = [reshape(dummyidx( 1:ui-1, uj+1:end), (ui-1)*(n-uj),1); ...
            reshape(dummyidx(ui+1:end, 1:uj-1), (n-ui)*(uj-1), 1)];
        dis(tableRowIdx, :, :) = sum( table(idx, :,:), 1);
    end
    
    conSum = squeeze(sum(table.*con, 1)) ./ 2;
    disSum = squeeze(sum(table.*dis, 1)) ./ 2;
    
    gamma = (conSum - disSum) ./ (conSum + disSum);
    
%     s2 = squeeze(sum(table .* ((con-dis).^2), 1)) - ((conSum-disSum).^2)./sampleSize;
%     se = 2./(conSum+disSum) * sqrt(s2);
%         
%     z = gamma./se; %inference via z-score
%     
%     if tail == 0
%         p = 2*(1 - normcdf(abs(z))); %p-value (two-sided)
%     else
%         p = 1 - normcdf(abs(z)); %p-value (one-sided)
%     end