function res = jiggle(data, scaleofvar, noise)
    %jiggle data: add tiny noise to the data 
    %for plotting
    %
    %scale: factor to multiply the variance of the data, def = 1e-2
    %
    %
    
    if nargin < 2
        scaleofvar = 1e-1;
    end
    if nargin < 3
        noise = [];
    end
    
    if ~isempty(noise)
        res = data + randn(size(data)) .* noise;
    else
        v = nanvar(data);
        res = data + randn(size(data)) .* repmat(v,size(data,1),1) * scaleofvar;
    end
end