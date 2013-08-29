function [pval res] = permTest_beta(nrun, rseed, tail, fhandle, data, varargin)
% for permutation test
% nrun: #runs
% rseed: random seed, for reproducing results
% fhandle: function handle, eg. @nanCorr
% data: data to be randomly permuted; should be a column vector
% varargin: all other inputs to the function
%   
% revised: from permTest.m, rewrite for large number of permutation
%   instead of looping over nrun, loop over data series (column)
% revised: 03-22-08, rseed=-1 & memsave=true, generate permindices on the
%   fly; rseed ~= -1, generate permindices first
% revised: 09-25-07, tail for one- or two-tail pval calculation
%   'tail1' or 'tail2'
% revised: 10-08-07, tailr, taill, right or left tail
% revised: 12-20-07, matrix results, data can be matrix
%

    t = cputime();
    [nsamples, ncol] = size(data);
        
    baserun = 1000; %at least run 1K
    if nrun < baserun
        error('permTest_beta is to run at least 1K.\n');
    end
    reprun = nrun/baserun;    
      
    if rseed ~= -1
        rand('state',rseed);
    end    
    
    signal = fhandle(data, varargin{:});
    %nseries: #series to permute (permute samples in each series)
    %ntest: #target to test on
    [nseries, ntest] = size(signal);
    pval = NaN(nseries, ntest);
    
    res = signal;
    count = ones(size(res));
    
    for repi = 1:reprun
        permindices = zeros(nsamples,baserun);
        for basei = 1:baserun
            permindices(:,basei) = randsample(nsamples, nsamples);
        end
        for i = 1:nseries %should be the same as ncol        
            x = data(:,i);
            tmp = fhandle(x(permindices), varargin{:}); %baserun * ntest
            resforcmp = repmat(res(i,:),baserun,1);
            switch tail
                case 'tail1'
                    for j = 1:ntest
                        if res(i,j) > 0
                            count(i,j) = count(i,j) + sum( tmp(:,j)>=res(i,j) ) ;
                        else
                            count(i,j) = count(i,j) + sum( tmp(:,j)<=res(i,j) );
                        end
                    end
                case 'tailr'
                    count(i,:) = count(i,:) + sum(tmp >= resforcmp);
                case 'taill'
                    count(i,:) = count(i,:) + sum(tmp <= resforcmp);
                otherwise %tail2
                    count(i,:) = count(i,:) + sum(abs(tmp) >= abs(resforcmp));
            end
        end
    end
                 
    pval = count / (nrun+1);
    fprintf('time: %f\n',cputime()-t);
end
