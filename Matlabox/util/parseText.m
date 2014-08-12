function res = parseText(fn, varargin)
    %fn: file name
    %varargin:
    %   'skip': number of lines to skip at the beginning of the file
    %       (comments etc)
    %   'colname': {true, false}; def=false;
    %   'rowname': {true, false}; def=false; if true, rownames are extracted
    %       from first column
    %   'ncol': number of columns if known
    %   'delimiter': default=tab
    %   'numeric': if the text body is numeric (except for the rowname, if
    %       exists); default=false;
    %   'nrowname': default = 0; number of columns in the beginning should
    %       be read as the "name" for rows
    %   'ncolname': default = 0; number of rows in the the beginning of the
    %       file should be read as the "name" for columns
    %   'numericcol': vector to specify which columns are numeric (indices
    %   are counted from the first column, including row names
    %   
    %
    %Note: to make sure it reads file correctly, put NaN for ''
    %
    %return:
    %   res.rowname, if rowname=true
    %   res.colname, if colname=true
    %   res.text
    %
    
    %default
    delimiter = '\t'; %9; %tab
    colname = false;
    rowname = false;
    skip = 0;
    ncol = 0;
    numeric = false;
    nrowname = 0;
    ncolname = 0;
    bufsize = 4095;
    numericcol = [];
    commentstyle = '#';
    
    passpara = {};
    i = 1;
    while i < length(varargin)
        switch lower(varargin{i})
            case {'colname','rowname','numeric'}
                eval(sprintf('%s = logical(%d);',varargin{i},varargin{i+1}));
            case {'skip','ncol', 'nrowname', 'ncolname', 'bufsize'}
                eval(sprintf('%s = %d;',varargin{i},varargin{i+1}));
            case {'numericcol'}
                eval(sprintf('%s = %s;', varargin{i}, mat2str(varargin{i+1})));
            case {'commentstyle'}
                eval(sprintf('%s = ''%s'';', varargin{i}, varargin{i+1}));
            case {'delimiter'}
                if strcmp(varargin{i+1},'\t')
                    delimiter = 9;
                else
                    delimiter = varargin{i+1};
                end
            otherwise
                passpara(end+1:end+2) = varargin(i:i+1);
                %error('Unknown option %s.\n',varargin{i});
        end
        i = i + 2;
    end
    
    if rowname && nrowname == 0
        nrowname = 1;
    end
    if colname && ncolname == 0
        ncolname = 1;
    end
    if nrowname >= 1 && ~rowname
        rowname = true;
    end
    if ncolname >= 1 && ~colname
        colname = true;
    end
    
    if ~isempty(numericcol)
        if size(numericcol,1) > size(numericcol, 2)
            numericcol = numericcol';
        end
    end
    
    fid = fopen(fn);    
    %skip lines if any
    i = 0;
    line = -1;
    while i < skip
        line = fgetl(fid);
        i = i + 1;
    end
    offset = 0;
    if ~ischar(line)
        fpos = ftell(fid);
        line = fgetl(fid);
        offset = fpos - ftell(fid);        
    end    
    while ~isempty( regexp(line, sprintf('^%s',commentstyle), 'once'))
        fpos = ftell(fid);
        line = fgetl(fid);
        offset = fpos - ftell(fid);        
    end
    fseek(fid, offset, 'cof');        
        
    if colname        
        for i = 1:ncolname
            hline = fgetl(fid);            
            a = textscan(hline, '%s', 'delimiter', delimiter);
            if i == 1
                res.colname = cell(length(a{1}), ncolname);
            end
            res.colname(:, i) = a{1};
        end
    end
    
    if numeric && isempty(numericcol)
        formattail = ' %f';     res.text = [];
    else
        formattail = ' %s';     res.text = {};
    end
    
    format = '';
    if isnumeric(delimiter), delimiter = '\t'; end
    
    if ncol == 0
        %read one data line
        fpos = ftell(fid);
        line = fgetl(fid);
        if ~ischar(line) %nothing to read
            return
        end
        offset = fpos - ftell(fid);
        fseek(fid, offset, 'cof');
        %decide the number of columns
        a = textscan(line, '%s', 'delimiter', delimiter, passpara{:});
        if colname
            ncol = size(res.colname,1);
            %read a data line to confirm the number of columns; this is
            %because sometimes the header does not have labels for the
            %first column
            if ncol == length(a{1}) - 1
                ncol = ncol + 1;
                res.colname = [{''}; res.colname];
            end
        else            
            ncol = length(a{1});
        end
    end
    for i = 1:nrowname
        format = [format ' %s'];
    end
    for i = 1:ncol-nrowname
        if ismember(i+nrowname, numericcol)
            format = [format ' %f'];
        else
            format = [format formattail];
        end
    end        
    
    text = textscan(fid, format, 'delimiter',delimiter, 'BufSize', bufsize, passpara{:});    
    fclose(fid);
    
    removecolnameidx = [];
    if rowname
        res.rowname = cell(length(text{1}), nrowname);
        for i = 1:nrowname
            res.rowname(:,i) = text{i};
        end
        text = text(nrowname+1:end);
        if colname
            res.rownamelabel = res.colname( 1:nrowname );
            removecolnameidx = [removecolnameidx 1:nrowname];
        end
    end
    
    nrows = length(text{1});
    if ~isempty(numericcol)        
        res.numtext = NaN(nrows, length(numericcol));
        for i = 1:length(numericcol)
            res.numtext(:,i) = text{ numericcol(i) - nrowname };            
        end        
        textcol = setdiff(1:length(text), numericcol-nrowname);        
        res.text = cell(nrows, length(textcol));
        for i = 1:length(textcol)
            res.text(:,i) = text{ textcol(i) };
        end
        if colname
            res.numcolname = res.colname(numericcol);
            removecolnameidx = [removecolnameidx numericcol];            
        end
    else
        if iscell(text{1})
            res.text = cell(nrows, length(text));
        else
            res.text = NaN(nrows, length(text));
        end
        for i = 1:length(text)
            res.text(:,i) = text{i};
        end
    end

    if colname
        res.colname( removecolnameidx, : ) = [];
    end
    
    