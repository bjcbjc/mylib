function gridforimage(varargin)


    x = xlim;
    y = ylim;
    
    xtick = x(1)+1 : (x(2)-1); 
    ytick = y(1)+1 : (y(2)-1);
    
    linex = [ repmat(xtick, 2, 1), repmat(x', 1, length(ytick))]; 
    liney = [ repmat(y', 1, length(xtick)), repmat(ytick, 2, 1) ];

    hold on;
    line( linex, liney, varargin{:});
    hold off;
    