function [nr, nc] = numSubplot(nplots)

b = sqrt(nplots);
l = floor(b);
u = ceil(b);
r = round(b);

if r == u
    nr = u;
    nc = u;
else
    nr = l;
    nc = u;
end
