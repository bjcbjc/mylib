function cvidx = cvtestidx(n, ncv)
    %return vector 1 x n: each element is the batch of test
    %eg: 5-cv: vector of n element with number 1-5 randomly distributed,
    %indicating which samples are test samples for each cv run
    cvidx = [repmat(1:ncv, 1, floor(n / ncv)), 1:mod(n, ncv)];
    cvidx = cvidx(randperm(n));
end