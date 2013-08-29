function [beta p] = fastcorr(x, y, n)
    %x, y: training entries and targets
    %y: #sample * #feature
    %x: #sample * #feature
    %
    %
 
    
    %this is original
%     xsqsum = sum(x.^2)'; %column            
%     xty = x'*y;
%     beta = (xty) ./ (xsqsum); %column
     
    %but if x and y are normalized, sum(x(:,1).^2) = n-1
    beta = x'*y ./ (n-1);

    if nargout > 1
        p = 2*tcdf(-abs(beta.*sqrt((n-2)./(1-beta.^2))),n-2);
%         p = ones(size(beta));
%         k = beta ~= 0;
%         p(k) = 2*tcdf(-abs(beta(k).*sqrt((n-2)./(1-beta(k).^2))),n-2);
    end
end