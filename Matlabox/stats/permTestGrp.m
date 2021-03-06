function [pval res] = permTestGrp(nrun, rseed, tail, fhandle, memsave, data, group, varargin)
% for permutation test
% nrun: #runs
% rseed: random seed, for reproducing results
% fhandle: function handle, eg. @nanCorr
% data: data to be randomly permuted; should be a column vector
% memsave: logical, if false, res returns all the results from permutation
%   runs; if true, res only return the real signal
% varargin: all other inputs to the function
%
% revised: 03-23-08, copied from permTest; add option for group (genotype)
%   so that permutation is restricted within each group
% revised: 03-22-08, rseed=-1 & memsave=true, generate permindices on the
%   fly; rseed ~= -1, generate permindices first
% revised: 03-20-08, add in memsave option
% revised: 09-25-07, tail for one- or two-tail pval calculation
%   'tail1' or 'tail2'
% revised: 10-08-07, tailr, taill, right or left tail
% revised: 12-20-07, matrix results, data can be matrix
    
    [nsamples, ncol] = size(data);    
    uv = unique(group);
    if length(uv) ~= 2, error('genotype error\n'); end
    
    if isscalar(rseed)
        if rseed ~= -1
            rand('state',rseed);
            permindices = zeros(nsamples,nrun);                            
            for ri = 1:nrun
                for geno = 1:length(uv)
                    k = find(group==uv(geno));
                    permindices(k,ri) = randsample(k,length(k));
                end                
            end
        else
            permindices = -1;
        end
    else        
        if size(rseed,2) ~= nrun || size(rseed,1) ~= nsamples
            error('rseed matrix size is not consistant, %dx%d\n',...
                size(rseed,1),size(rseed,2));
        end
        permindices = rseed;
    end
    
    signal = fhandle(data, varargin{:});
    %nseries: #series to permute (permute samples in each series)
    %ntest: #target to test on
    [nseries, ntest] = size(signal);
    pval = NaN(nseries, ntest);
    if ~memsave        
        res = NaN(nseries, ntest, nrun+1);
        res(:,:,1) = signal;
        for ri = 1:nrun
            res(:,:,ri+1) = fhandle(data(permindices(:,ri),:), varargin{:});
        end        
        switch tail
            case 'tail1'        
                for i = 1:nseries
                    for j = 1:ntest
                        if res(i,j,1) > 0
                            pval(i,j) = length(find(squeeze(res(i,j,:))>=res(i,j,1))) / (nrun+1);
                        else
                            pval(i,j) = length(find(squeeze(res(i,j,:))<=res(i,j,1))) / (nrun+1);
                        end
                    end
                end
            case 'tailR'
                for i = 1:nseries
                    for j = 1:ntest
                        pval(i,j) = length(find(squeeze(res(i,j,:))>=res(i,j,1))) / (nrun+1);
                    end
                end
            case 'tailL'
                for i = 1:nseries
                    for j = 1:ntest
                        pval(i,j) = length(find(squeeze(res(i,j,:))<=res(i,j,1))) / (nrun+1);
                    end
                end
            case 'tail2'
                for i = 1:nseries
                    for j = 1:ntest
                        pval(i,j) = length(find(abs(squeeze(res(i,j,:)))>=abs(res(i,j,1)))) / (nrun+1);
                    end
                end
            otherwise
                error('Unknown tail %s\n',tail);
        end
    else
        res = signal;
        count = ones(size(res));
        
        for ri = 1:nrun
            if isscalar(permindices)                
                permi = zeros(nsamples,1);
                for geno = 1:length(uv)
                    k = find(group==uv(geno));
                    permi(k) = randsample(k,length(k));
                end    
                tmp = fhandle(data(permi,:), varargin{:});
            else
                tmp = fhandle(data(permindices(:,ri),:), varargin{:});
            end
            switch tail
                case 'tail1'
                    for i = 1:nseries
                        for j = 1:ntest                            
                            if res(i,j) > 0
                                count(i,j) = count(i,j) + ( tmp(i,j)>=res(i,j) ) ;
                            else
                                count(i,j) = count(i,j) + ( tmp(i,j)<=res(i,j) );                                
                            end
                        end
                    end
                case 'tailR'
                    count = count + (tmp >= res);
                case 'tailL'
                    count = count + (tmp <= res);                    
                case 'tail2'
                    count = count + ( abs(tmp) >= abs(res) );                    
                otherwise
                    error('Unknown tail %s\n',tail);
            end
        end
        pval = count / (nrun+1);
    end
    %fprintf('time: %f\n',cputime()-t);
end
