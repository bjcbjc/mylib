function [rowidx, colidx] = hierOrder(mtx, distance, xy)

if nargin < 2,
    distance = 'correlation';
end

if nargin < 3,
    xy = [true, true];
end

rowidx = [];
colidx = [];
if xy(1)
    d = pdist(mtx, distance);
    tree = linkage(d, 'average');
    rowidx = optimalleaforder(tree, d);    
end

if xy(2)
    d = pdist(mtx', distance);
    tree = linkage(d, 'average');
    colidx = optimalleaforder(tree, d);    
end