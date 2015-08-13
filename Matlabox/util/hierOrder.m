function [rowidx, colidx] = hierOrder(mtx, distance, xy, optimal)

if nargin < 2,
    distance = 'correlation';
end

if nargin < 3,
    xy = [true, true];
end

if nargin < 4,
    optimal = true;
end

rowidx = [];
colidx = [];
if xy(1)
    d = pdist(mtx, distance);
    tree = linkage(d, 'average');
    if optimal
        rowidx = optimalleaforder(tree, d);    
    else
        [~, ~, rowidx] = dendrogram(tree, 0);
    end
end

if xy(2)
    d = pdist(mtx', distance);
    tree = linkage(d, 'average');
    if optimal
        colidx = optimalleaforder(tree, d);    
    else
        [~, ~, colidx] = dendrogram(tree, 0);
%         colidx = tree(:,1:2);
%         colidx = colidx(colidx <= size(mtx, 2));
    end
end