function [pval res] = permTest(nrun, rseed, tail, fhandle, data, varargin)

% for permutation test
% nrun: #runs
% rseed: random seed, for reproducing results
% fhandle: function handle, eg. @nanCorr
% data: data to be randomly permuted; should be a column vector
% varargin: all other inputs to the function
%   

% revised: 09-25-07, tail for one- or two-tail pval calculation
%   'tail1' or 'tail2'
% revised: 10-08-07, tailr, taill, right or left tail

    %t = cputime();
    [ndata, ncol] = size(data);
    if ncol > 1
        error('data is not a column vector\n');
    end
    
    if isscalar(rseed)
        rand('state',rseed);
        permindices = zeros(ndata,nrun);
        for i = 1:nrun
            permindices(:,i) = randsample(ndata,ndata);
        end
    else        
        if size(rseed,2) ~= nrun || size(rseed,1) ~= ndata
            error('rseed matrix size is not consistant, %dx%d\n',...
                size(rseed,1),size(rseed,2));
        end
        permindices = rseed;
    end
    
    %res = zeros(nrun+1,1);    
    %res(1) = fhandle(data, varargin{:});
    %res(2:end) = fhandle(data(permindices), varargin{:});
    signal = fhandle(data, varargin{:});
    ntest = size(signal,2);
    res = zeros(nrun+1,ntest);
    res(1,:) = signal;
    res(2:end,:) = fhandle(data(permindices), varargin{:});    
    pval = zeros(1,ntest);
    
    if strcmp(tail, 'tail1')
        for i = 1: ntest
            if res(1,i) > 0
                pval(i) = length(find(res(:,i)>=res(1,i))) / (nrun+1);
            else
                pval(i) = length(find(res(:,i)<=res(1,i))) / (nrun+1);
            end
        end
    elseif strcmp(tail, 'tailr')
        for i = 1:ntest
            pval(i) = length(find(res(:,i)>=res(1,i))) / (nrun+1);
        end
    elseif strcmp(tail, 'taill')
        for i = 1:ntest
            pval(i) = length(find(res(:,i)<=res(1,i))) / (nrun+1);
        end
    else %tail2
        for i = 1:ntest
            pval(i) = length(find(abs(res(:,i))>=abs(res(1,i)))) / (nrun+1);
        end
    end
    %fprintf('time: %f\n',cputime()-t);
end