function generateGeneSet(fn, genes, conds, membership, fileformat)
    % output to fn
    % genes are the first column, ORF
    % conds are the first row
    % membership is #conds cells, each cell contains a vector
    % of indices corresponding to genes
    % fileformat: {'tab','xml'}
    %
    %
        
    if strcmp(fileformat, 'tab')
        gset_values = zeros(length(genes),length(conds));

        fid = fopen(fn,'w');
        fprintf(fid,'Gene');
        for i = 1:length(conds)
            if ~isempty(strfind(conds{i},'&'))
                fprintf(fid,'\t%s',strrep(conds{i},'&','+'));
            else
                fprintf(fid,'\t%s',conds{i});
            end
            %set bit value by the way
            gset_values(membership{i}, i) = 1;
        end
        fprintf(fid,'\n');
        for i = 1:length(genes)
            if sum(gset_values(i,:)) == 0 %skip genes that are not assigned
                continue
            end
            fprintf(fid,'%s',genes{i});
            fprintf(fid,'%s\n',num2str(gset_values(i,:),'\t%d'));
        end
        fclose(fid);
    elseif strcmp(fileformat, 'xml')
        fid = fopen(fn, 'w');
        fprintf(fid, '<Root>\n');
        for i = 1:length(conds)
            fprintf(fid, '<Module Name="%s">\n',strrep(conds{i},'&','+'));
            fprintf(fid, '<Set>\n');
            for j = 1:length(membership{i})
                fprintf(fid, '%s\t', genes{membership{i}(j)});
            end
            fprintf(fid, '\n</Set>\n');
            fprintf(fid, '</Module>\n');
        end
        fprintf(fid, '</Root>\n');
        fclose(fid);
    else
        fprintf('unknown format %s\n',fileformat);
    end
end

