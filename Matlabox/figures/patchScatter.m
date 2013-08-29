function dothandle = patchScatter(x, y, varargin)

    para.facealpha = 0.3;
    para.facecolor = 'b';
    para.edgecolor = 'none';
    para.size = 10;
    para.nvertex = 5;
    para.xsize = [];
    para.ysize = [];
    
    para = assignpara(para, varargin{:});
    
    if numel(para.facecolor) > 1
        if numel(para.facecolor) ~= numel(x)
            [m, n] = size(para.facecolor);
            if n == 3
                para.cdata = reshape(para.facecolor, 1, m, n);
            elseif m == 3
                para.cdata = reshape(para.facecolor', 1, n, m);
            else
                para.cdata = para.facecolor;
            end
        else
            para.cdata = para.facecolor;
        end
        para.facecolor = 'flat';
    end
    
    figpos = get(gcf, 'position');
    axspos = get(gca, 'position');
    
    if isempty(para.xsize) || isempty(para.ysize)
        xx = xlim;
        yy = ylim;
        if xx(1) > min(x) || xx(2) < max(x)
            xlim([min(x) max(x)]);
        end
        if yy(1) > min(y) || yy(2) < max(y)
            ylim([min(y) max(y)]);
        end
        xsize = range(xlim) / figpos(3) / axspos(3) * para.size / 2;
        ysize = range(ylim) / figpos(4) / axspos(4) * para.size / 2;
    else
        xsize = para.xsize;
        ysize = para.ysize;
    end    
        
    t = (0:pi/para.nvertex:2*pi)';
    
    para = rmfield(para, {'size', 'nvertex','xsize','ysize'});
    
    x = x(:)';
    y = y(:)';
    xdata = bsxfun(@plus, xsize.*sin(t), x);
    ydata = bsxfun(@plus, ysize.*cos(t), y);

    paras = paraPairArray(para);
        
    hold on;
    dothandle = patch( 'xdata', xdata, 'ydata', ydata);
    set(dothandle, paras{:} );
    hold off;

    if para.facealpha < 1 && ~strcmpi(get(gcf, 'renderer'), 'opengl')
        set(gcf, 'renderer', 'opengl');
    end
end