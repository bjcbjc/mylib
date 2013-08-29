function acc = weightedAccuracy(label, pred, weight, classLabels)

    %label: n x nY
    %pred: n x nY
    %weight: nClass x nY, or 1 x nY cell
    %classLabels: vector or cell

    if nargin < 4
        classLabels = [];
    end
    
    [~, nY] = size(label);
        
    if isempty(classLabels)
        classLabels = cell(1, nY);
        for i = 1:nY
            classLabels{i} = mynanunique(label(:,i));            
        end
        nClass = cellfun(@length, classLabels);
    else
        if iscell(classLabels)
            nClass = cellfun(@length, classLabels);
        else
            nClass = length(classLabels);
        end
    end
    
    if iscell(weight)
        nClassWeight = cellfun(@length, weight);        
    else
        nClassWeight = size(weight, 1);        
    end
    if any(nClass > nClassWeight)
        fprintf('not enough weight for #class!\n');
        acc = [];
        return
    end
    
    acc = zeros(1, nY);
    if iscell(classLabels)
        for i = 1:nY
            subacc = zeros(nClass(i) , 1);
            for k = 1:nClass(i)
                subacc(k) = sum(label(:,i) == classLabels{i}(k) & pred(:,i) == classLabels{i}(k)) ./ ...
                    sum(label(:,i) == classLabels{i}(k));
            end
            acc(i) = sum(weight{i}.*subacc) ./ sum(weight{i});
        end
    else
        subacc = zeros(nClass, 1);
        for i = 1:nY
            subacc(:) = 0;
            for k = 1:nClass
                subacc(k) = sum(label(:,i) == classLabels(k) & pred(:,i) == classLabels(k)) ./ ...
                    sum(label(:,i) == classLabels(k));
            end
            acc(i) = sum(weight(:,i) .* subacc) ./ sum(weight(:,i));
        end
    end
    
end