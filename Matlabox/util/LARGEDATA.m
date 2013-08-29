classdef LARGEDATA < handle

    
    methods (Static)
        function s = nanstd(data, flag, dim)
            idx = isnan(data);
            n = size(data, dim) - sum(idx, dim);
            if flag == 1
                df = n;
            else
                df = max(n-1, 1);
            end
            data(idx) = 0;
            m0 = sum(data, dim) ./ n;
            res = bsxfun(@minus, data, m0);
            res(idx) = 0;
            s = sqrt( sum(res .^ 2, dim) ./df );
        end
        
        function s = nanvar(data, flag, dim)
            idx = isnan(data);
            n = size(data, dim) - sum(idx, dim);
            if flag == 1
                df = n;
            else
                df = max(n-1, 1);
            end
            data(idx) = 0;
            m0 = sum(data, dim) ./ n;
            res = bsxfun(@minus, data, m0);
            res(idx) = 0;
            s = sum(res .^ 2, dim) ./df ;
        end
        
        function r = binMtx(datafn, rowIdx, colIdx, mtxSize)
            %
            
            n = length(rowIdx);
            r = NaN(n, 1);
            idx = sub2ind(mtxSize, rowIdx, colIdx);
            
            segments = find( diff(idx) > 1);            
            segments = [0; segments; n];
            f = fopen(datafn, 'r');
            for i = 1:length(segments)-1
                nToRead = segments(i+1) - segments(i);
                start = idx( segments(i) + 1 ) * 4;
                fseek(f, start, 'bof');
                data = fread(f, nToRead, 'single');
                
            for k = 1:length(usubcolidx)
                byteidx = (sub2ind( [blocksize(uniblockidx(j)) ...
                    blocksize(uniblockidx(i)) ], m, usubcolidx(k) ) - 1) .* 4;
                fseek(f, byteidx, 'bof');
                columndata = fread(f, M-m+1, 'single');
                r(r_rowidx, r_colidx( ucoli == k ) ) = columndata(submtx_rowidx-m+1);
            end
            fclose(f);
            
            
            
        end
        
        function r = retrivePairCorr(datahead, indices, nfeature, nblock, ifmtx)
            % extract values from a HUGE corr matrix which has been broken down to
            % blocks (by blockcorr)
            %
            % datahead: directory/file_header for each block of the corr matrix
            % indices: list of indices indexing the original feature list
            % nfeature: number of features (the complete list)
            % nblock: number of blocks the corr matrix is broken into
            % ifmtx: return matrix or vector
            %
            % The function will extract the submatrix (#indices x #indices) from
            % the complete corr matrix. If ifmtx is false, the function vectorizes
            % the lower triangle submatrix
            %
            
            if nargin < 5
                ifmtx = false;
            end
            
            nidx = length(indices);
            r = NaN(nidx);
            
            % decide which blocks need to be read
            [blockidx, withinidx, blocksize] = fidx2blockidx(indices, nfeature, nblock);
            [uniblockidx, ~, uidx] = unique(blockidx);
            for i = 1:length(uniblockidx)
                r_colidx = find(uidx == i);
                submtx_colidx = withinidx( r_colidx );
                [usubcolidx, ~, ucoli] = unique(submtx_colidx);
                for j = i:length(uniblockidx)
                    fn = [datahead sprintf('.%02d.%02d.bin', uniblockidx(j), uniblockidx(i))];
                    r_rowidx = find(uidx == j);
                    submtx_rowidx = withinidx( r_rowidx );
                    m = min( submtx_rowidx );
                    M = max( submtx_rowidx );
                    f = fopen(fn, 'r');
                    for k = 1:length(usubcolidx)
                        byteidx = (sub2ind( [blocksize(uniblockidx(j)) ...
                            blocksize(uniblockidx(i)) ], m, usubcolidx(k) ) - 1) .* 4;
                        fseek(f, byteidx, 'bof');
                        columndata = fread(f, M-m+1, 'single');
                        r(r_rowidx, r_colidx( ucoli == k ) ) = columndata(submtx_rowidx-m+1);
                    end
                    fclose(f);
                end
            end            
            
            if ifmtx
                i = isnan(r);
                rt = r';
                r(i) = rt(i);
            else
                i = tril(ones(nidx));
                r = r( i == 1);
            end
        end
        
        function [blockidx withinidx blocksize] = fidx2blockidx(fidx, nfeature, nblock)
            k = floor( nfeature / nblock );
            alpha = mod( nfeature, nblock);
            blocksize = k .* ones(1, nblock);
            blocksize(1:alpha) = blocksize(1:alpha) + 1;
            csize = cumsum(blocksize);
            blockidx = arrayfun(@(x) find( csize >= x, 1, 'first'), fidx);
            csize = [0 csize];
            withinidx = arrayfun(@(x) fidx(x) - csize( blockidx(x) ), 1:length(fidx))';
        end
        
    end
    
end