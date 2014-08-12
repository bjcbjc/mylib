function prepareGSEAFile(expression, phenotype, outPath, outHead)
    %expression: struct, with .data, .sample, .name
    %phenotyp: struct, with .data, .sample, .name
    %
    %

    %align sample first
    [~, si, sj] = intersect(expression.sample, phenotype.sample);

    expression.sample = expression.sample(si);
    expression.data = expression.data(:, si);
    phenotype.sample = phenotype.sample(sj);
    phenotype.data = phenotype.data(:, sj);

    if ~iscell(phenotype.name)
        phenotype.name = {phenotype.name};
    end
    
    %output expression
    %
    fn = sprintf('%s/%s.exp.txt', outPath, outHead);
    dlmwriteplus(fn, expression.data, ...
        [expression.name, repmat({'NA'}, length(expression.name), 1)], ...
        [{'gene', 'description'}, expression.sample(:)']);

    if iscell(phenotype.data)
        n = length(phenotype.sample);
        for i = 1:length(phenotype.name)
            fn = sprintf('%s/%s.phenotype.%s.categorical.cls', outPath, outHead, phenotype.name{i});
            f = fopen(fn, 'w');                    
            uLabel = unique(phenotype.data(i,:), 'stable');
            fprintf(f, '%d %d 1\n', n, length(uLabel));
            fprintf(f, '#%s\n', sprintf(' %s', uLabel{:}));
            fprintf(f, '%s\n', sprintf('%s ', phenotype.data{i,:}));
            fclose(f);
        end        
    else
        fn = sprintf('%s/%s.phenotype.numeric.cls', outPath, outHead);
        f = fopen(fn, 'w');
        fprintf(f, '#numeric\n');
        for i = 1:length(phenotype.name)
            fprintf(f, '#%s\n', phenotype.name{i});
            fprintf(f, '%s\n', sprintf('%g ', phenotype.data(i,:)));
        end
        fclose(f);
    end
    
    
    