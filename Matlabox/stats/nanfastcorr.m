function [beta] = nanfastcorr(x, y, upattern, ui)
    %x, y: training entries and targets
    %y: #sample * 1, a column vector
    %x: #sample * #feature
    %
    %
 
    %x is normalized
    yvi = ~isnan(y');
    x = x';
    
    
    
    if nargin < 4
        xvi = ~isnan(x);        
        [upattern, ~, ui] = unique(xvi, 'rows');
        upattern = bsxfun(@and, upattern, yvi);
        clear xvi 
    end
    
    beta = NaN(length(ui), 1);
    x = zeromean_univar_normalization(x, 2);
    for i = 1:size(upattern, 1)
        %vi = yvi & upattern(i,:); 
        %if sum(vi) == 0, continue; end
        if sum( upattern(i,:) ) == 0, continue; end
        nrmy = zeromean_univar_normalization(y( upattern(i,:) ), 1);
        
        beta(ui==i) = x(ui==i, upattern(i,:))*nrmy ./ (sum(upattern(i,:))-1);
    end
    

    
end