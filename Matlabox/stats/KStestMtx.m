function [p] = KStestMtx(data, label)
%An interface to take matrix or vector data to do Kolmogorov-Smirnov test
%Each column in data is taken as a series of data and categorized by
%label, which can either be a single vector or a matrix. Number of rows
%in label must be the same as that in data

    [nsample, nseries] = size(data);
    [tmp, ncat] = size(label);
    if nsample == 1 && nseries > 1
        data = data';
        [nsample, nseries] = size(data);
    end
    
    if tmp == 1 && ncat > 1
        label = label';
        [tmp, ncat] = size(label);
    end
    t = cputime();
    if nsample ~= tmp
        error('numbers of samples are different in data and label.\n');
    end
    
    p = NaN(nseries, ncat);
    for i = 1:nseries
        for j = 1:ncat
            unilabel = unique(label(:,j));
            if length(unilabel) ~= 2
                error('more than 2 types of labels are input.\n');
            end
            x1 = data(find(label(:,j)==unilabel(1)),i);
            x2 = data(find(label(:,j)==unilabel(2)),i);           
            [h, p(i,j)] = kstest2(x1,x2);
        end
    end        
    fprintf('time:%f\n',cputime()-t);
end
