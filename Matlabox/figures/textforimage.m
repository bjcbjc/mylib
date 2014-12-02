function h = textforimage(textdata, varargin)
    para.numformat = '%0.2f';
    para.fontsize = 10;
    para.color1 = [0, 0, 0];
    para.color2 = [0.7, 0.7, 0.7];
    para.idxForColor2 = [];
    para.pctColor2 = 20;
    para.noNaN = true;
    para.color2Dir = 'both';
       
    para = assignpara(para, varargin{:});
    [x, y] = meshgrid(1:size(textdata,2), 1:size(textdata,1));
    if isnumeric(textdata)
        mtxstring = strtrim(cellstr(num2str(textdata(:), '%0.2f')));
    else
        mtxstring = textdata;
    end
    if para.noNaN
        mtxstring = strrep(mtxstring, 'NaN', '');
    end
    h = text(x(:), y(:), mtxstring(:), 'horizontalalignment', 'center', ...
        'fontsize', para.fontsize, 'color', para.color1);
        
    if ~isempty( para.idxForColor2)
        set(h(para.idxForColor2), 'color', para.color2);
    elseif isnumeric(textdata)
        n = size(colormap(), 1);
        a = get(colorbar('peer', gca), 'ytick');
        a = a(1):(a(end)-a(1))/(n-1):a(end);
        thres1 = a(max(1, floor( n * para.pctColor2/100)));
        thres2 = a(min(n, floor( n * (100-para.pctColor2)/100)));
        if strcmp(para.color2Dir, 'both')
            set(h( textdata <= thres1 | textdata >= thres2 ), 'color', para.color2);
        elseif strcmp(para.color2Dir, 'high')
            set(h( textdata >= thres2 ), 'color', para.color2);
        else strcmp(para.color2Dir, 'low')
            set(h( textdata <= thres1 ), 'color', para.color2);
        end
    end
