function gloc = index2gloc(index)

    maxbp = 300e6; %max bp ~250e6, chr1
    
    gloc = [ floor(index/maxbp)+1, mod(index, maxbp) ];
    
end