function dec = binvec2dec(bindata)

    [n, m] = size(bindata);
    dec = zeros(n, 1);
    p = (m-1):-1:0;
    for i = 1:m
        dec = dec + (2^p(i)) .* bindata(:,i);
    end
