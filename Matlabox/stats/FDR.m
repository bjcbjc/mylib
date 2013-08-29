function [pc, q] = FDR(qc, p)
% qc: q-value control
% p: p-values, a vector
% Benjamini and Hochberg procedure
%
% return:
% q: q-value
% pc: cutoff for p-val corresponding to q
%

    if length(p) == 1 || length(qc) > 1
        error('usage: FDR(qc, p)\n');        
    end

    valididx = ~isnan(p);
       
    [r c] = size(p(valididx));
    [sp pi] = sort(p(valididx));
    
    if r>c
        k = (1:r)'; m = r;
    else
        k = 1:c;    m = c;
    end
    
    %pc = sp( max( find( sp <= qc/m*k ) ) );
    pc = sp( find( sp <= qc/m*k, 1, 'last') );
    if isempty(pc)
        pc = 0;
    end
    if nargout > 1
        valididx = find(valididx);
        for i = m-1:-1:1
            sp(i) = min(sp(i+1), sp(i)*m/i);
        end
        q = NaN(size(p));
        q(valididx(pi)) = sp;
    end
end