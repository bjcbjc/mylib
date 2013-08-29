function res = strkey(varargin)
    % str_or_num: string, cell array of string, or number(s)
    % varargin can be {'query', 'remove', 'savelater', 'savenow','add','saveto','lookup','dbname'}
    % 
    % strkey(str), strkey(strs), strkey(strs, 'query') return unique numeric keys
    %   for input strings; new keys are generated if strings are not found
    %   in current look-up table and 'add' flag is on; otherwise, return 0
    %   for those not in the table
    % strkey(num) returns string or cell array of strings represented by
    %   the key 'num'; empty string is returned for those not in the table
    % strkey(str, 'remove') removes string(s) from lookup table
    % strkey(..., 'savelater') don't save the table if new keys are added;
    %   but users should save the table later if the table is to be used in
    %   the future. This is to save I/O time when consenquent addition of
    %   new keys are expected
    % strkey('savenow') saves the current table; default behavior is to
    %   save the table if any new keys are added
    % strkey(..., 'saveto', filename)
    % strkey(..., 'lookup', filename)
    % 'saveto', 'lookup': specify the filename of .mat; error if both flags
    %   are specified and are not the same
    % strkey('dbname'): return the current mat file name used
    %
    
    persistent STRKEY_HASH;
    persistent matname; 
        
    flag = {'remove', 'savelater', 'savenow', 'add', 'saveto', 'lookup','dbname'};
    flagidx = cellfun(@ischar, varargin);
    flagvalue = ismember(flag, varargin(flagidx));
    flagidx = find(flagidx);
    
    for i = 1:length(flag)
        eval(sprintf('%s=flagvalue(i);', flag{i}));
    end
    
    ifsave = true;    
    if ~savenow && savelater
        ifsave = false;    
    end
    if saveto && lookup
        savetoidx = flagidx( strcmp(varargin(flagidx), 'saveto') ) + 1;
        lookupidx = flagidx( strcmp(varargin(flagidx), 'lookup') ) + 1;
        if ~strcmp(varargin{savetoidx}, varargin{lookupidx})
            error('saveto %s and lookup %s do not agree', varargin{savetoidx}, varargin{lookupidx});
        end
        matname = varargin{savetoidx};
    elseif saveto
        savetoidx = flagidx( strcmp(varargin(flagidx), 'saveto') ) + 1;
        matname = varargin{savetoidx};
    elseif lookup
        lookupidx = flagidx( strcmp(varargin(flagidx), 'lookup')) + 1;
        matname = varargin{lookupidx};
    end
    
    if isempty(matname)
        matname = 'STRKEY_HASH.mat';
    end
    
    dbdir = strrep(mfilename('fullpath'), mfilename, '');
    dbfilename = [dbdir matname];
        
    if isempty(STRKEY_HASH)
        if exist(dbfilename, 'file')
            STRKEY_HASH = loadStructData(dbfilename);
        else
            STRKEY_HASH.id = 0;
            STRKEY_HASH.key = {''};
            STRKEY_HASH.needsave = true;
            STRKEY_HASH.matname = matname;
        end
    end
    
    if ~strcmpi(STRKEY_HASH.matname, matname) 
        if lookup %change lookup table
            if exist(dbname, 'file')
                STRKEY_HASH = loadStructData(dbfilename);
            else
                STRKEY_HASH.id = 0;
                STRKEY_HASH.key = {''};
                STRKEY_HASH.needsave = true;
                STRKEY_HASH.matname = matname;
            end
        else %change saved filename
            STRKEY_HASH.needsave = true;
        end
    end

    if nargin == 0
        if isempty(STRKEY_HASH)
            fprintf('No look-up table is available\n');
        else
            fprintf('Current STRKEY_HASH: %d keys\n', length(STRKEY_HASH.key));
        end
        res = []; 
        return
    end
    
    if dbname
        res = STRKEY_HASH.matname;
    elseif remove
        if ischar(varargin{1})
            varargin{1} = {varargin{1}};
        end
        tf = ismember(STRKEY_HASH.key, varargin{1});
        STRKEY_HASH.key(tf) = [];
        STRKEY_HASH.id(tf) = [];
        STRKEY_HASH.needsave = true;
    elseif ~savenow
        if ~isempty(STRKEY_HASH.key)
            lastid = STRKEY_HASH.id(end);
        else
            lastid = 0;
        end
        str_or_num = varargin{1};        
        if ischar(str_or_num)
            str_or_num = {str_or_num};
        end
        
        n = length(STRKEY_HASH.key);
        if iscellstr(str_or_num)             
             [tf, res] = ismember(str_or_num, STRKEY_HASH.key);
             if any(~tf(:)) && add %new key
                 [ukey, ~, uidx] = unique(str_or_num(~tf));
                 naddkey = length(ukey);
                 addidx = n+1:n+naddkey;                 
                 STRKEY_HASH.key(addidx,1) = ukey;
                 STRKEY_HASH.id(addidx,1) = lastid+1:lastid+naddkey;
                 STRKEY_HASH.needsave = true;
                 res(~tf) = addidx(uidx);
                 tf(~tf) = true;
             end
             res(tf) = STRKEY_HASH.id(res(tf));
        else
            [tf, idx] = ismember(str_or_num, STRKEY_HASH.id);
            if any(~tf(:))
                res = cell(size(str_or_num));
                res(:) = {''};
                res(tf) = STRKEY_HASH.key(idx(tf));
            else
                res = STRKEY_HASH.key(idx);
            end
        end        
    end
    
    
    if STRKEY_HASH.needsave && ifsave
        STRKEY_HASH.needsave = false;
        save(dbfilename, 'STRKEY_HASH');
    end
    
    