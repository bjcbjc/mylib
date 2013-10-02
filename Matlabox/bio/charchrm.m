function chrm = charchrm(chrm, padchr, mitoformat)
    if nargin < 2
        padchr = true;
    end
    if nargin < 3
        mitoformat = 'M';
    end
    
    chrm = numarray2strarray(chrm);
    chrm(strcmp(chrm, '23')) = {'X'};
    chrm(strcmp(chrm, '24')) = {'Y'};
    chrm(strcmp(chrm, '25')) = {mitoformat};
    
    if padchr
        chrm = strcat('chr', chrm);
    end
end