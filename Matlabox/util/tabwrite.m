function tabwrite(fname, cellarray, writeAsOneString)
    %write to file with fname, tab-delimited
    %
    %
    if nargin < 3, writeAsOneString = false; end
    [r, c] = size(cellarray);
    
    if r/c >= 100 && ~writeAsOneString
        writeAsOneString = true;
        fprintf('many lines with few columns, use writeAsOneString option\n');
    elseif c >= 10 && r >= 1000 && writeAsOneString
        fprintf('many columns, disable writeAsOneString option\n');
        writeAsOneString = false;
    end
    
    if writeAsOneString         
        formatstr = repmat({'%s'}, 1, c);
        formatstr( cellfun(@isnumeric, cellarray(1,:)) ) = {'%g'};
        formatstr = strjoin(formatstr, '\\t');
        cellarray = cellarray';
        fid = fopen(fname, 'w');
        fprintf(fid, [formatstr '\n'], cellarray{:});
        fclose(fid);
    else
        fid = fopen(fname, 'w');
        try
            if iscellstr(cellarray)
                for i = 1:r
                    fprintf(fid, '%s\n', strjoin(cellarray(i,:), '\t'));
                end    
            else
                for i = 1:r
                    formatstr = '';
                    valstr = '';
                    for j = 1:c
                        if isnumeric(cellarray{i,j})
                            formatstr = strcat(formatstr, '%g\t');
                        else
                            formatstr = strcat(formatstr, '%s\t');
                        end
                        valstr = strcat(valstr, sprintf(',cellarray{%d,%d}',i,j));
                    end
                    eval(sprintf('fprintf(fid, ''%s\\n'' %s);', formatstr, valstr));
                end
            end
            fclose(fid);
        catch
            fclose(fid);
            rethrow(lasterror);
            error('error; no file is written');
        end
    end
end