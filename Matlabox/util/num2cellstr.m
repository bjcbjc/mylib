function cstr = num2cellstr(numarray)

[m, n] = size(numarray);
cstr = reshape( strtrim( cellstr( num2str( numarray(:) ) ) ), m, n);
