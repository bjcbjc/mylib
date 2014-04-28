function idx = jobDivider(jobidx, numTotalJobs, n)
    %jobidx = index of job
    %n: number of tasks
    k = floor(n / numTotalJobs);
    alpha = mod(n, numTotalJobs);
    if jobidx <= alpha
        idx = (jobidx-1) * (k+1) + 1;
        idx = idx : idx + k;
    else
        idx = alpha * (k+1) + (jobidx-1-alpha) * k + 1;
        idx = idx : idx + (k-1);
    end
end