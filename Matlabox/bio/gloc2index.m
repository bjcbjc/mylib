function index = gloc2index(gloc, gloc2)
    if nargin < 2
        gloc2 = [];
    end
    
    maxbp = 300e6; %max bp ~250e6, chr1
    
    if ~isempty(gloc2)
        if ischar(gloc)
            gloc = numericchrm(gloc);
        end
        index = (gloc-1) * maxbp + gloc2;
    else
        index = (gloc(:,1)-1) * maxbp + gloc(:,2);
    end