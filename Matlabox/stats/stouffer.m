function [combinedP, combinedZ] = stouffer(z, w)
    % z: #z-stat x #meta-test
    % weights for each z-stat: #z-stat x #meta-test
    %
    % combine Z scores with Stouffer's method:  z = sum(zi *wi) / sqrt(sum
    %   wi^2)
    % assume the z-stats for each meta-test are independent
    
    if nargin < 2, w = ones(size(Z)); end
    
    combinedZ = sum( z .* w, 1) ./ sqrt(sum(w.^2, 1));
    combinedP = 2*normcdf(abs(combinedZ), 'upper');
    