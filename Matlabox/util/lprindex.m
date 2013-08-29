function index = lprindex(i, j, n, inv, triangle)
    % For large data, symmetric matrix can be saved with only upper or
    % lower triangle data. The symmetric values in the matrix imply the
    % pair-wise relationships between two data vectors is the same regardless of
    % the order of the two data vectors. For example, corr(i,j) =
    % corr(j,i). Therefore, we can save only data(i,j) with i<j, with
    % linearized index. Assume the data is stored in the order of (i,j)
    % where i<j, from j = 1:n, for each i (i=1~n-1), this function is to
    % return the linearized index to access the data stored.
    %
    
    if nargin < 4, inv = false; end
    if nargin < 5, triangle = 'lower'; end
    
    if ~inv
        assert(length(i)==length(j), '#i ~= #j');

        if strcmpi(triangle, 'lower')
            i1 = min(i, j);
            j1 = max(i, j);
            
            assert(all(i1<=n), 'i1 <= n does not hold');
            assert(all(j1<=n), 'j1 <= n does not hold');
            
            index = ( (2*n - i1) .* (i1-1) )./2 + (j1-i1);
            index( i1 == j1 ) = 0;
        else
            i1 = min(i, j);
            j1 = max(i, j);
            
            assert(all(i1<=n), 'i1 <= n does not hold');
            assert(all(j1<=n), 'j1 <= n does not hold');
            
            index = (j1-1) .* (j1-2) ./2 + i1;
            index( i1 == j1 ) = 0;
            
        end
    else
        %only i is valid; i should be the indices to the linearized vector;
        %here we inverse i into (i,j) to the original symmetric matrix
        %n is also need; but j is ignored
        if size(i, 2) > size(i, 1)
            i = i';
        end
        k = length(i);
        index = zeros(k, 2);
        
        if strcmpi(triangle, 'lower')
            ncumsum = cumsum( (n-1 : -1 : 1)' );
            index(:, 2) = arrayfun(@(x) find(ncumsum >= x, 1, 'first'), i);
            ncumsum = [0; ncumsum];
            index(:, 1) = i - ncumsum(index(:, 2)) + index(:, 2);
        else
            ncumsum = cumsum( (1:(n-1))' );
            index(:, 2) = arrayfun(@(x) find(ncumsum >= x, 1, 'first'), i);
            ncumsum = [0; ncumsum];
            index(:, 1) = i - ncumsum(index(:, 2));
            index(:, 2) = index(:, 2) + 1;
        end
        
%         for j = 1:k
%             cur = i(j);
%             rowcount = n-1;
%             index(j, 1) = 1;
%             while cur > rowcount
%                 cur = cur - rowcount;
%                 index(j,1) = index(j,1) + 1;
%                 rowcount = rowcount - 1;
%             end
%             index(j, 2) = index(j, 1) + cur;
%         end
    end