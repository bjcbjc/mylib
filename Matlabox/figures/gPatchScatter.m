function h = gPatchScatter(x, y, g, varargin)

    para.renderOrder = '';
    para.facealpha = 0.3;
    para.facecolor = [];
    para.edgecolor = 'none';
    para.size = 10;
    para.nvertex = 5;
    
    defcolor = {'r','g', 'b', 'y', 'm', 'c', 'o', 'p', 'n'};
    defn = length(defcolor);
    
    para = assignpara(para, varargin{:});
    
    figpos = get(gcf, 'position');
    axspos = get(gca, 'position');
    xx = xlim;
    yy = ylim;
    if xx(1) > min(x) || xx(2) < max(x)
        xlim([min(x) max(x)]);
    end
    if yy(1) > min(y) || yy(2) < max(y)
        ylim([min(y) max(y)]);
    end
    para.xsize = range(xlim) / figpos(3) / axspos(3) * para.size / 2;
    para.ysize = range(ylim) / figpos(4) / axspos(4) * para.size / 2;
    
    
    [ugroup, ~, groupidx] = unique(g);
        
    ngroup = length(ugroup);
    if isempty(para.facecolor)
        para.facecolor = ColorDict.translate( defcolor(mod((1:ngroup)-1, defn)+1));
    elseif isscalar(para.facecolor)
        m = size(para.facecolor, 1);
        if m == 3 && ngroup ~= 3 %rgb
            para.facecolor = para.facecolor';
        end
    end
    
    sidx = 1:ngroup;
    if ~isempty(para.renderOrder)
        if iscellstr(para.renderOrder) && length(para.renderOrder) == length(ugroup)
            [~, orderNum] = ismember(ugroup, para.renderOrder);
            [~, sidx] = sort(orderNum);
        end        
    end

    allFaceColor = para.facecolor;
    
    h = zeros(ngroup, 1);
    for i = 1:ngroup
        idx = groupidx == sidx(i);
        para.facecolor = allFaceColor(i,:);
        parastr = paraPairArray(para);
        h(i) = patchScatter(x(idx), y(idx), 'PASSPARA', parastr{:});
        if iscellstr(ugroup)
            set(h(i), 'displayname', ugroup{sidx(i)});
        else
            set(h(i), 'displayname', num2str(ugroup(sidx(i))));
        end
    end
    legend('show');
    
        
    