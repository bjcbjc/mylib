function chrm = numericchrm(chrm)

    chrm = strrep(chrm, 'chr', '');
    chrm = strrep(chrm, 'X', '23');
    chrm = strrep(chrm, 'Y', '24');
    chrm = strrep(chrm, 'MT', '25');
    chrm = strrep(chrm, 'M', '25');

    chrm = str2double_fast(chrm);
end