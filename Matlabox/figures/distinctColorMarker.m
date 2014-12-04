function handles = distinctColorMarker(handles, varargin)
    %  assign color and marker attributes for each handle so that they can
    %  be identified; for example, it can be used gscatter plot when the
    %  number of groups is large
    
    para.line = {'-', '--', ':', '-.'};
    para.marker = {'.', 'sq', '^', 'o', 'x', '*'};
    para.color = [0.9, 0.1, 0.1; ... %red
        0.1, 0.8, 0.1; ... %green
        0.1, 0.1, 0.9; ... %blue        
        1, 0, 1; ... %magenta
        0.1, 0.8, 0.8; ... %cyan
        0.98, 0.67, 0.15; ... %orange        
        0.3, 0.3, 0.3; ...%black
        0.7, 0.7, 0.7]; %grey
    para.markersize = [15, 8, 8, 8, 8];
    para.maxcolor = size(para.color, 1);
    para.useLine = false;
    para.loopColorFirst = true;
    
    para = assignpara(para, varargin{:});
    if iscellstr(para.color)
        para.color = ColorDict.translate(para.color);
    end
    
    ns = length(para.markersize);
    nm = length(para.marker);
    nl = length(para.line);
    if ns < nm
        para.markersize(ns+1:nm) = para.markersize(ns);
    end

    if para.maxcolor < size(para.color, 1)
        para.color = para.color(1:para.maxcolor, :);
    end
    
    nc = size(para.color, 1);    
    nh = length(handles);
    if para.useLine
        if nc*nl < nh
            error('not enough color * line');
        end
        if para.loopColorFirst
            [cidx, midx] = ind2sub([nc, nl], 1:nh);
        else
            [midx, cidx] = ind2sub([nl, nc], 1:nh);
        end
        for i = 1:nh
            if isprop(handles(i), 'color')
                set(handles(i), 'color', para.color(cidx(i), :));
            end
            if isprop(handles(i), 'linestyle')
                set(handles(i), 'linestyle', para.line{midx(i)});
            end
        end
    else
        if nc*nm < nh
            error('not enough color * marker');
        end
        if para.loopColorFirst
            [cidx, midx] = ind2sub([nc, nm], 1:nh);
        else
            [midx, cidx] = ind2sub([nm, nc], 1:nh);
        end
        for i = 1:nh
            if isprop(handles(i), 'color')
                set(handles(i), 'color', para.color(cidx(i), :));
            end
            if isprop(handles(i), 'marker')
                set(handles(i), 'marker', para.marker{midx(i)}, 'markersize', para.markersize(midx(i)));
            end
        end
    end
    
        