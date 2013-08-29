function [pval res] = permTest(nrun, rseed, fhandle, data, varargin)

% for permutation test
% nrun: #runs
% rseed: random seed. for reproducing results
% fhandle: function handle, eg. @nanCorr
% data: data to be randomly permuted
% varargin: all other inputs to the function
%     

    t = cputime();
    rand('state',rseed);
    res = zeros(nrun+1,1);
    res(1) = fhandle(data, varargin{:});
    ndata = length(data);
    for i = 1:nrun
        res(i+1) = fhandle(data(randsample(ndata,ndata)), varargin{:});
    end
    
    if res(1) > 0
        pval = length(find(res>=res(1))) / (nrun+1);
    else
        pval = length(find(res<=res(1))) / (nrun+1);
    end
    fprintf('time: %f\n',cputime()-t);
end