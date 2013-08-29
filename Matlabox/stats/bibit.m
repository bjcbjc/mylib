function [clusterRows, clusterCols] = bibit(X, minr, minc)

nr = sum(X,1);
nc = sum(X,2);
vc = nr >= minr;
vr = nc >= minc;

X = X(vr, vc);
[n1, n2] = size(X);
transposed = false;
if n2 < n1
    X = X';
    [n1, n2] = size(X);
    tmp = minr;
    minr = minc;
    minc = tmp;
    tmp = vr;
    vr = vc;
    vc = tmp;
    transposed = true;
end
%n2 > n1, to reduce n2


wordsize = 32;
bins = [1:wordsize:n2, n2+1];
numBins = length(bins)-1;
code = zeros(n1, numBins); %n1 x reduced_dim
for i = 1:numBins
    code(:,i) = bin2dec(num2str(X(:, bins(i):bins(i+1)-1),'%1d'));
end

lastpattern = 0;
cRows = false(n1, 1000);
cCols = false(n2, 1000);
rowMembership = false(n1, 1);
for i = 1:n1-1
    seeds = bsxfun(@bitand, code(i,:), code(i+1:n1,:));
    seeds = unique(seeds, 'rows'); %uidx to recover membership of rows
        
    nCmember = sum(reshape(sum( ...
        dec2bin(seeds)=='1',2),size(seeds)),2); %#seeds x 1, how many cols in each pattern    
    seeds(nCmember < minc, :) = [];
                
    for si = 1:size(seeds,1)
        rowMembership(:) = all( bsxfun(@eq, ...
            bsxfun(@bitand, seeds(si,:), code), seeds(si,:)),2);
        
        if sum(rowMembership) < minr, continue; end
                
        colMembership = (dec2bin(seeds(si,:)) == '1');
        part1 = reshape( ...
            [ zeros(numBins-1, wordsize-size(colMembership,2)), ...
            colMembership(1:end-1,:) ]', (numBins-1)*wordsize, 1);
        nAddZero = bins(end)-bins(end-1)-size(colMembership,2);
        if nAddZero >= 0
            colMembership = [ part1; ...
                zeros(nAddZero,1); colMembership(end,:)' ];
        else
            colMembership = [ part1; ...
                colMembership(end,end-(bins(end)-bins(end-1))+1:end)' ];
        end
        
        %is a subset?
        if ~any(sum( bsxfun(@and, cRows(:,1:lastpattern), rowMembership),1) ... 
                == sum(rowMembership) & ...
                sum( bsxfun(@and, cCols(:,1:lastpattern), colMembership), 1) ...
                == sum(colMembership))
            lastpattern = lastpattern + 1;
            cRows(:,lastpattern) = rowMembership;
            cCols(:,lastpattern) = colMembership;
        end
        
        if lastpattern >= size(cRows,2)
            cRows = [cRows false(n1,1000)];
            cCols = [cCols false(n2,1000)];
        end
        
    end
end

remove = ~any(cRows,1);
cRows(:, remove) = [];
cCols(:, remove) = [];
nclusters = sum(~remove);

if ~transposed
    clusterRows = false(length(vr), nclusters);
    clusterCols = false(length(vc), nclusters);
    clusterRows(vr,:) = cRows;
    clusterCols(vc,:) = cCols;
else
    clusterRows = false(length(vc), nclusters);
    clusterCols = false(length(vr), nclusters);
    clusterRows(vc,:) = cCols;
    clusterCols(vr,:) = cRows;
end


