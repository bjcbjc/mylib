function [rho p] = spearman_perm(X, Y, numperm, randomseed, permstep)
    %X: #sample x #variable
    %Y: #sample x #variable
    %numperm
    %
    
    
    if nargin < 3
        numperm = 1000;
    end
    if nargin < 4
        randomseed = 100;
    end
    if nargin < 5
        permstep = 10000;
    end
    
    if any(isnan(X(:))) || any(isnan(Y(:)))
        error('cannot have NaN');        
    end
    
    %set up blocks of permutations    
    if numperm < permstep
        perms = numperm;
    else
        perms = [repmat(permstep, 1, floor(numperm/permstep)), mod(numperm,permstep)];
    end
        
    perms(perms==0) = [];
    
    [n, nx] = size(X);
    ny = size(Y,2);
    
    if ny < nx
        tmp = Y;
        Y = X;
        X = tmp;
        swapped = true;
        clear tmp
        nx = size(X,2);
        ny = size(Y,2);
    else
        swapped = false;
    end
    
    rho = NaN(nx, ny);
    p = zeros(nx, ny);
    
    n3const = (n+1)*n*(n-1) ./ 3;
    
    %temp use
    Xadj = NaN(1, nx);
    Yadj = NaN(1, ny);
    meanD = NaN(nx, ny);
    stdD = NaN(nx, ny);
    
    for i = 1:ny
        [Y(:,i), Yadj(i)] = tiedrank(Y(:,i), 0);
    end
    for i = 1:nx
        [X(:,i) Xadj(i)] = tiedrank(X(:,i), 0);
        D = sum(bsxfun(@minus, X(:,i), Y).^ 2, 1); %1 x nY
        meanD(i,:) = (n3const - (Xadj(i) + Yadj)./ 3) ./ 2; %1 x nY
        stdD(i,:) = sqrt( (n3const./2 - Xadj(i)./3) .* (n3const./2 - Yadj./3)./(n-1) ); %1 x nY
        rho(i, :) = (meanD(i,:) - D) ./ (sqrt(n-1)*stdD(i,:)); %1 x nY
    end

    setrandomseed(randomseed);
    %permutation, no need to calculate the rank again, but permute
    %the rank; only D will change since "adj" will not be changed
    
    randidx = NaN(n, perms(1)); %n x perms(pi)
    for pi = 1:length(perms)
        if pi == length(perms)
            randidx(:, perms(pi)+1:end) = [];
        end
        for ii = 1:perms(pi)
            randidx(:, ii) = randperm(n);
        end
        for xi = 1:nx
            x = X(:,xi); %n x 1
            for yi = 1:ny
                D = sum( bsxfun(@minus, x( randidx), Y(:, yi)).^2, 1); %1 x perms(pi)
                permrho = (meanD(xi,yi) - D) ./ (sqrt(n-1)*stdD(xi,yi)); %1 x perms(pi)
                p(xi, yi) = p(xi, yi) + sum(permrho >= rho(xi,yi));
            end
        end
    end
    
    p = p ./ numperm;
    
    if swapped
        rho = rho';
        p = p';
    end
    
end
                
%     [xrank, xadj] = tiedrank(xi,0);
%     [yrank, yadj] = tiedrank(yj,0);
%     D = sum((xrank - yrank).^2);
%     meanD = (n3const - (xadj+yadj)./3) ./ 2;
%     stdD = sqrt((n3const./2 - xadj./3)*(n3const./2 - yadj./3)./(n-1));
%     coef(i,j) = (meanD - D) ./ (sqrt(n-1)*stdD);