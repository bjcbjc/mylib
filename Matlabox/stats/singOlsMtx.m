function [beta] = singOlsMtx(x, y)
    %x, y: training entries and targets
    %y: #sample * 1, a column vector
    %x: #sample * #feature
    %
    %
    %It can be used to calculate correlation coefficient if x and y are
    %standardized.
    %

    
    %[nTrainSamples, nFeatures] = size(x);
            
    xsqsum = sum(x.^2)'; %column            
    xty = x'*y;
    %w = sum((my-w0).*x,2) ./ (sum(x.^2,2)); %a column of w
    beta = (xty) ./ (xsqsum); %column
     
    % #sample * #feature
    %residual = repmat(y,1,nFeatures) - repmat(beta',nTrainSamples,1).*x;
   
    %faster for large x
    %residual = bsxfun(@minus, y, bsxfun(@times, beta', x));

    
end