function datastruct = affy2gene(fn, dictionary, fnout)
    %fn is data after mas5;
    %first line is annotation
    %then two columns as [probe id, intensity]
    %
    %dictionary is a file with two colums: [probe id, gene id]
    %
    %
    
    data = parseText(fn, 'skip', 1, 'ncol', 2, 'nrowname', 1, 'ncolname', 0, 'numeric', true);
    dict = parseText(dictionary, 'ncol', 2, 'nrowname', 1, 'ncolname', 0, 'numeric', true);
    
    [~, idx] = ismember(data.rowname, dict.rowname);
    data.rowname(idx == 0) = [];
    data.text(idx == 0, :) = [];
    data.rowname = dict.text(idx(idx~=0));
    
    if ~isempty(fnout)
        dlmwriteplus(fnout, data.text, data.rowname, {});
    end    
    datastruct.gene = data.rowname;
    datastruct.value = data.text;

    
    

