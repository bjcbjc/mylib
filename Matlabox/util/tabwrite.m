function tabwrite(fname, cellarray)
    %write to file with fname, tab-delimited
    %
    %
    
    [r c] = size(cellarray);
    
    fid = fopen(fname, 'w');
    try
        if iscellstr(cellarray)
            for i = 1:r
                fprintf(fid, '%s\n', strjoin(cellarray(i,:), '\t'));
            end
%             formatstr = '%s';
%             for i = 2:c
%                 formatstr = strcat(formatstr, '\t%s');
%             end
%             formatstr = [formatstr '\n'];
%             for i = 1:r
%                 valstr = '';
%                 for j = 1:c
%                     valstr = strcat(valstr, sprintf(',cellarray{%d,%d}',i,j));
%                 end
%                 eval(sprintf('fprintf(fid, ''%s'' %s);', formatstr, valstr));
%             end
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