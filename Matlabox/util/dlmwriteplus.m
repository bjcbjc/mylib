function dlmwriteplus(filename, mtx, rowname, colname, delim, varargin)

    if ~isempty(rowname) && isrow(rowname) && size(mtx,1) > 1
        rowname = rowname';
    end
    if ~isempty(colname) && iscolumn(colname) && size(mtx,2) > 1
        colname = colname';
    end
    
    [row_nr, row_nc] = size(rowname);
    [col_nr, col_nc] = size(colname);
    [mtx_nr, mtx_nc] = size(mtx);
    
    if ~isempty(rowname)
        assert(row_nr==mtx_nr, '#rowname ~= size(mtx,1)');
    end
    if ~isempty(colname)
        assert(col_nc==mtx_nc || col_nc == mtx_nc+row_nc, ...
            '#colname ~= size(mtx,2) or size(mtx,2)+size(rowname,2)');
    end
    if nargin < 5
        delim = sprintf('\t');
    end
    
    ntotalcol = mtx_nc + row_nc;
    ntotalrow = mtx_nr + col_nr;
    
    %padding
    if ntotalcol > col_nc
        colname = [repmat({''}, col_nr, ntotalcol - col_nc), colname];
    end
    
    f = fopen(filename, 'w');

    if ~isempty(colname)        
        for ri = 1:col_nr
            fprintf(f, '%s\n', strjoin(colname(ri,:), delim));
        end
    end
    
    if isempty(rowname)
        fclose(f);
        if ~isempty(colname)
            dlmwrite(filename, mtx, '-append', 'delimiter', delim, varargin{:});
        else
            dlmwrite(filename, mtx, 'delimiter', delim, varargin{:});
        end
    else
%         for ri = 1:row_nr
%             fprintf(f, '%s\n', strjoin(rowname(ri,:), delim));
%         end
%         dlmwrite(filename, mtx, '-append', 'roffset', col_nr, 'coffset', row_nc, 'delimiter', delim);
        tmpfn = sprintf('tmp%04d.txt', randi(1000));
        dlmwrite(tmpfn, mtx, 'delimiter', delim, varargin{:});
        ftmp = fopen(tmpfn, 'r');
        
        li = 1;
        line = fgets(ftmp);
        while ischar(line)
            fprintf(f, '%s%s%s', strjoin(rowname(li,:),delim), delim, line);
            line = fgets(ftmp);
            li = li + 1;
        end
    
        fclose(ftmp);
        fclose(f);
        
        system(sprintf('rm -f %s',tmpfn));
    end
                    
    
    