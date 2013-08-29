function res = breakbiclik(E, coli, colj, minr, minc) 
    %input E: #row > #col
    %break down E for each coli:colj, taking the submatrix that has all
    %ones in each column in coli:colj
    %    
    %input to biclique should have larger number of rows than columns
    %
    
    runcol = colj - coli + 1;
    res.col = (coli:colj)';
    res.minr = minr;
    res.minc = minc;
    res.bic = cell(runcol, 2);    
    
    for runi = 1:runcol
        subE = E( E(:,coli+runi-1) == 1, :); %cut dow #rows by restriction in a column 
        [sr sc] = size(subE);
        if sr < sc %need to transpose
            bic = biclique(subE', 'minu', minc, 'minv', minr);
            res.bic{runi, 1} = bic{2};
            res.bic{runi, 2} = bic{1};
        else
            bic = biclique(subE, 'minu', minr, 'minv', minc);
            res.bic{runi, 1} = bic{1};
            res.bic{runi, 2} = bic{2};
        end
    end
    
    