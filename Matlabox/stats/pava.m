function Yest = pava(Y, x)
    %Y: n x #test
    %x: 1 x n or n x 1
    [x, indx] = sort(x, 'ascend');
    [~, indxinv] = sort( indx );
    
    Y = Y(indx, :);
    n = length(indx);
    nY = size(Y, 2);
    
    lvlsets = (1:n)';
    xunique = unique(x);
    
    for i = 1:length(xunique)
       valssamex = x == xunique(i);
       Y(valssamex, :) = repmat(mean( Y(valssamex, :), 1), sum(valssamex), 1);           
       lvlsets(valssamex) = min( lvlsets(valssamex) );
    end
        
    lvlsets = repmat(lvlsets, 1, nY);
    while true
      viol = ( Y(1:(n-1), :) - Y(2:n, :) > 0 );
      yi = find( any(viol, 1) )';
      if isempty(yi)
          break; 
      end
      
      [~, i] = sort(viol(:, yi), 1, 'descend');
      i = i(1,:)';
      lvl1 = lvlsets( sub2ind( [n nY], i, yi ) )';
      lvl2 = lvlsets( sub2ind( [n nY], i+1, yi ) )';
      ilvl = bsxfun(@eq, lvlsets(:, yi), lvl1) | bsxfun(@eq, lvlsets(:, yi), lvl2);
      yhat = sum( Y(:, yi) .* ilvl, 1) ./ sum( ilvl, 1);
      yhat = repmat(yhat, n, 1);
      yhat = yhat(ilvl);
      [ri, ci] = find(ilvl);
      indx = sub2ind([n nY], ri, yi(ci) );
      Y( indx ) = yhat;
      yhat = repmat(lvl1, n, 1);      
      lvlsets( indx ) = yhat( sub2ind( [n length(yi)], ri, ci) );      
    end
    
    Yest = Y(indxinv, :);    
end

