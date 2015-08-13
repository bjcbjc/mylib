function [H, AX, BigAX, P, PAx, textH] = plotmatrix_symm(X, label, varargin)

    [H, AX, BigAX, P, PAx] = plotmatrix(X);
    set(P, 'facecolor', [0.7, 0.7, 0.9], 'edgecolor', 'w');
    
    m = min(X(:));
    M = max(X(:));
    d = (M - m) * 0.05;
    m = m - d;
    M = M + d;
    
    set(AX, 'xlim', [m, M], 'ylim', [m, M]);
    n = size(X, 2);
    textH = NaN(1, n);
    for i = 1:n
        axes(PAx(i));
        set(AX(i,i), 'xticklabel', get(AX(i,i), 'xtick'), ...
            'yticklabel', get(AX(i,i), 'ytick'));
        x = xlim;
        y = ylim;
        textH(i) = text(x(1)+diff(x)*0.5, y(1)+diff(y)*0.5, label{i}, ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle', ...
            varargin{:});
        delete(AX(i+1:end, i));
    end
    
    r = corr(X, 'rows','pairwise');
    rho = corr(X, 'rows','pairwise', 'type', 'spearman');
    for i = 1:n-1
        for j = i+1:n
            axes(AX(i,j));
            x = xlim;
            y = ylim;
            text(x(1)+diff(x)*0.1, y(2)-diff(y)*0.1, ...
                sprintf('r = %0.2f\nrho = %0.2f', r(i,j), rho(i,j)), ...
                'fontsize', 12, ...
                'verticalalignment', 'top');
        end
    end