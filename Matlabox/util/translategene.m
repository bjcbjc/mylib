function res = translategene(db, querylist, from, to)
    if nargin < 4 
        to = 'gene';
    end
    if nargin < 3
        from = 'orf';
    end
    
    if isempty(db)
        db = load('SGD');
        db = db.SGD;
    end
    
    assert(isfield(db, from), '%s is not in db', from);
    assert(isfield(db, to), '%s is not in db', to);
    
    if strcmpi(from, 'gene')
        fprintf('gene name might not be unique\n');
    end
    
    if ~iscellstr(querylist)
        querylist = {querylist};
    end


    res = querylist;
    
    [uquery ui uj] = unique(querylist);
    [tmp i j] = intersect(db.(from), uquery);
    uquery(j) = db.(to)(i);
    for i = 1:length(uquery)
        res(uj==i) = repmat(uquery(i), sum(uj==i), 1);
    end
    