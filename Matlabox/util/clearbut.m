function clearbut(varargin)
    %clear every variable but the input
    %take * wildcard but not ?
    
    vnames = evalin('base', 'whos');
    vnames = {vnames.name};
    n = length(vnames);
    
    keep = false(n, 1);

    haswild = find(cellfun(@(x) ismember('*',x), varargin));
    nowild = setdiff(1:length(varargin), haswild);
    
    for i = 1:length(nowild)
        keepi = strmatch(varargin{nowild(i)}, vnames, 'exact');
        if ~isempty(keepi)
            keep(keepi) = true;
        end
    end
    
    for i = 1:length(haswild)
        keepi = ~cellfun(@isempty, strfind(vnames, strrep(varargin{haswild(i)}, '*', '')));
        keep(keepi) = true;
    end
    

    toclear = vnames(~keep);
    clearstr = '';
    for i = 1:length(toclear)
        clearstr = sprintf('%s %s',clearstr, toclear{i});
    end
    if ~isempty(clearstr)
        clearstr = sprintf('clear %s',clearstr);
        evalin('base', clearstr);
    end
    
    