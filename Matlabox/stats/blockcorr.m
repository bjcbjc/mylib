function [r p] = blockcorr(data, blockperdim, blockidx1, blockidx2, usefastcorr, varargin)
    %data is large; calculate correlation for each block in a job
    %data should be #observation x #feature
    % blockperdim: number of blocks to divide #features
    % blockidx1, blockidx2: indices of blocks to the huge corr matrix (not
    %   calculated)
    % usefastcorr: use fastcorr function, but will not have p-value for
    %   spearman; if false, call corr(); if any NaN in the data, use corr()
    % varargin: pass to corr()
    %
    %The entire matrix should need #blockperdim * (#blockperdim+1) / 2 jobs
    %to finish.
    %
    
    if nargin < 5 || any(isnan(data(:)))
        usefastcorr = false;
    end
            
    [n p] = size(data);
    
%     if blockidx1 < blockidx2
%         tmp = blockidx2;
%         blockidx2 = blockidx1;
%         blockidx1 = tmp;
%     end
%     
    rowidx = jobDivider(blockidx1, blockperdim, p);
    colidx = jobDivider(blockidx2, blockperdim, p);
    if ~usefastcorr        
        if nargout > 1
            [r p] = corr(data(:, rowidx), data(:, colidx), varargin{:});
        else
            r = corr(data(:, rowidx), data(:, colidx), varargin{:});
        end
    else
        i = find(strcmpi(varargin, 'type'));
        ifspearman = strcmpi(varargin{i+1}, 'spearman');
        
        if ifspearman
            X = zeromean_univar_normalization( tiedrank( data(:, rowidx), 0), 1);
            if blockidx1 == blockidx2
                Y = X;
            else
                Y = zeromean_univar_normalization( tiedrank( data(:, colidx), 0), 1);
            end
            clear data
            r = fastcorr(X, Y, n);
        else
            X = zeromean_univar_normalization( data(:, rowidx), 1);
            if blockidx1 == blockidx2
                Y = X;
            else
                Y = zeromean_univar_normalization( data(:, colidx), 1);
            end
            clear data
            if nargout > 1
                [r p] = fastcorr(X, Y, n);
            else
                r = fastcorr(X, Y, n);
            end
        end
        
    end
    
    
    
    
    