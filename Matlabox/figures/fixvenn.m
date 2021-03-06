function [handles, setnum] = fixvenn(setdata, varargin)
    %setdata: #data x #set, should be logical (or will be converted to
    %   logical)
    %setdata: cell arrays; each cell represents data of one set;
    %   overlaps between sets will be computed and return in setnum
    %   structure
    %
    %options:
    %   'label': cell array, names of sets for display
    %   'color': rgb array; each row for each set
    %   'alpha': control tranparency of the diagram, def = 0.4
    %   'numfontsize': font size for the numbers of overlap between sets,
    %       def = 16
    %	'labelfontsize': font size for the labels of sets, def = 16
    %
    
    [~, i] = ismember('precalculated', varargin(1:2:end));
    [~, j] = ismember('setflag', varargin(1:2:end));
    i = max(2*i - 1, 0);
    j = max(2*j - 1, 0);
    if i ~= 0 
        precalculated  = varargin{i+1};
    else
        precalculated = false;
        if j ~= 0
            varargin(j:j+1) = [];
        end
    end
    if ~precalculated
        if iscell(setdata)
            u = setdata{1};
            for i = 2:length(setdata)
                u = union(u, setdata{i});
            end
            newsetdata = false(length(u), length(setdata));
            for i = 1:length(setdata);
                newsetdata(:,i) = ismember(u, setdata{i});            
            end
            setdata = newsetdata;
            clear newsetdata
        end
        if ~islogical(setdata)
            setdata = setdata ~= 0;
        end
        nset = size(setdata,2);
    else
        if ~isnumeric(setdata)
            error('precalculated data must be numeric');
        end
        setdata = setdata(:)';
        if length(setdata) == 3
            nset = 2;
        elseif length(setdata) == 7
            nset = 3;
        elseif length(setdata) == 15
            nset = 4;
        else 
            nset = -1;
        end
    end    
    
    
    %para = assignpara(para, varargin{:});
    switch nset
        case 2
            [handles, setnum] = venn2(setdata, varargin{:});
        case 3
            [handles, setnum] = venn3(setdata, varargin{:});
        case 4
            [handles, setnum] = venn4(setdata, varargin{:});
        otherwise
            fprintf('%d sets are not supported\n', nset);
    end
    
end

function [handles, setnum] = venn2(setdata, varargin)
    para.label = repmat({''}, 1, 2);
    para.color = [0.7, 0.3, 0.4;  0.4, 0.3, 0.8];
    para.alpha = 0.4;
    para.numfontsize = 16;
    para.labelfontsize = 16;
    para.precalculated = false;
    para.setflag = [0,1; 1,0; 1,1];
    para.setflagOrder = para.setflag;
    
    para = assignpara(para, varargin{:});
    
    if nnz(para.setflag ~= para.setflagOrder) && para.precalculated
        d1 = numarray2strarray(para.setflag);
        d2 = numarray2strarray(para.setflagOrder);
        k1 = arrayfun(@(x) strjoin(d1(x,:), '_'), 1:size(para.setflag, 1), 'unif', 0);
        k2 = arrayfun(@(x) strjoin(d2(x,:), '_'), 1:size(para.setflag, 1), 'unif', 0);
        [~, si] = ismember(k2, k1);
        para.setflag = para.setflag(si,:);
        setdata = setdata(si);
    end
    
    t = (0:2*pi/39:2*pi);    
    offset = [0, 0; -1, 0];
    labeloffset = [0.1, -0.1; -0.15, -0.1];
    labelanchorpoint = [16, 25];
    numlabelpos = [0.5,0; -1.5,0; -0.5,0];    
    
    x = 1;
    y = 1;
    xd = x*sin(t);
    yd = y*cos(t);
    
    cla
    hold on;
    %loc = get(gcf, 'position');
    %set(gcf, 'position', [loc(1), loc(2), 750, 650], 'color', 'w');
    handles.patch = NaN(2, 1);
    handles.label = NaN(2, 1);
    handles.numlabel = NaN(3, 1);
    
    for i = 1:2
        handles.patch(3-i) = patch(xd+offset(i,1), yd+offset(i,2), 'k');
        set(handles.patch(3-i), 'facecolor', para.color(i,:), 'facealpha', para.alpha, 'edgecolor', 'none');
        
        anchor = [xd(labelanchorpoint(i))+offset(i,1)+labeloffset(i,1), ...
            yd(labelanchorpoint(i))+offset(i,2)+labeloffset(i,2)];
        handles.label(3-i) = text( anchor(1), anchor(2), ...
            para.label{3-i}, 'fontsize', para.labelfontsize);        
        twidth = get(handles.label(3-i), 'extent');        
        if i == 2
            set(handles.label(3-i), 'position', [anchor(1) - twidth(3), ...
                anchor(2), 0]);        
        end
    end
    
    %fill in numbers of intersections
    setnum.set = para.setflag;
    setnum.num = zeros(size(para.setflag, 1), 1);    
    for i = 1:size(para.setflag,1)
        if ~para.precalculated
            setnum.num(i) = sum(all(bsxfun(@eq, setdata, para.setflag(i,:)), 2));
        else
            setnum.num(i) = setdata(i);
        end
        handles.numlabel(i) = text( numlabelpos(i,1), numlabelpos(i,2), ...
            num2str(setnum.num(i)), 'fontsize', para.numfontsize);
        twidth = get(handles.numlabel(i), 'extent');
        set(handles.numlabel(i), 'position', [numlabelpos(i,1) - twidth(3)/2, ...
            numlabelpos(i,2), 0]);
    end
    
    xlim([-2 1]);
    ylim([-1.5 1.5]);
    axis square 
    axis off
    hold off
end

function [handles, setnum] = venn3(setdata, varargin)
    para.label = repmat({''}, 1, 3);
    para.color = [0.7, 0.3, 0.4; 0.3, 0.7, 0.4; 0.4, 0.3, 0.8];
    para.alpha = 0.4;
    para.numfontsize = 16;
    para.labelfontsize = 16;
    para.precalculated = false;
    para.setflag = [0,0,1; 0,1,0; 1,0,0; 0,1,1; 1,0,1; 1,1,0; 1,1,1];
    para.setflagOrder = para.setflag;
    
    para = assignpara(para, varargin{:});
    
    if nnz(para.setflag ~= para.setflagOrder) && para.precalculated
        d1 = numarray2strarray(para.setflag);
        d2 = numarray2strarray(para.setflagOrder);
        k1 = arrayfun(@(x) strjoin(d1(x,:), '_'), 1:size(para.setflag, 1), 'unif', 0);
        k2 = arrayfun(@(x) strjoin(d2(x,:), '_'), 1:size(para.setflag, 1), 'unif', 0);
        [~, si] = ismember(k2, k1);
        para.setflag = para.setflag(si,:);
        setdata = setdata(si);
    end
    
    t = (0:2*pi/39:2*pi);    
    offset = [0, 0; -0.5, 1; -1, 0];
    labeloffset = [0.1, -0.1; 0, 0.1; -0.3, -0.1];
    labelanchorpoint = [16, 1, 25];
    numlabelpos = [0.20,-0.27; -0.76,1.36; -1.70,-0.27; -0.16,0.65; ...
        -0.76,-0.44; -1.36,0.65; -0.76,0.26];    
    
    x = 1;
    y = 1;
    xd = x*sin(t);
    yd = y*cos(t);
    
    cla
    hold on;
    %loc = get(gcf, 'position');
    %set(gcf, 'position', [loc(1), loc(2), 750, 650], 'color', 'w');
    handles.patch = NaN(3, 1);
    handles.label = NaN(3, 1);
    handles.numlabel = NaN(7, 1);
    
    for i = 1:3
        handles.patch(4-i) = patch(xd+offset(i,1), yd+offset(i,2), 'k');
        set(handles.patch(4-i), 'facecolor', para.color(i,:), 'facealpha', para.alpha, 'edgecolor', 'none');
        
        anchor = [xd(labelanchorpoint(i))+offset(i,1)+labeloffset(i,1), ...
            yd(labelanchorpoint(i))+offset(i,2)+labeloffset(i,2)];
        handles.label(4-i) = text( anchor(1), anchor(2), ...
            para.label{4-i}, 'fontsize', para.labelfontsize);        
        twidth = get(handles.label(4-i), 'extent');        
        if i == 3
            set(handles.label(4-i), 'position', [anchor(1) - twidth(3), ...
                anchor(2), 0]);
        elseif i == 2
            set(handles.label(4-i), 'position', [anchor(1) - twidth(3)/2, ...
                anchor(2), 0]);
        end
    end
    
    %fill in numbers of intersections
    setnum.set = para.setflag;
    setnum.num = zeros(size(para.setflag, 1), 1);    
    for i = 1:size(para.setflag,1)
        if ~para.precalculated
            setnum.num(i) = sum(all(bsxfun(@eq, setdata, para.setflag(i,:)), 2));
        else
            setnum.num(i) = setdata(i);
        end        
        handles.numlabel(i) = text( numlabelpos(i,1), numlabelpos(i,2), ...
            num2str(setnum.num(i)), 'fontsize', para.numfontsize);
    end
    
    xlim([-2.5 1.5]);
    ylim([-1.5 2.5]);
    axis square 
    axis off
    hold off
end

function [handles, setnum] = venn4(setdata, varargin)
    para.label = repmat({''}, 1, 4);
    para.color = [1, 0.5, 0.5; 0.5, 0.9, 0.5; 0.7, 0.6, 0.3; 0.5, 0.5, 1];
    para.alpha = 0.4;
    para.numfontsize = 16;
    para.labelfontsize = 16;
    para.precalculated = false;
    para.setflag = [0,0,0,1; 0,0,1,0; 0,0,1,1; 0,1,0,0; 0,1,0,1; 0,1,1,0; ...
        0,1,1,1; 1,0,0,0; 1,0,0,1; 1,0,1,0; 1,0,1,1; 1,1,0,0; 1,1,0,1; ...
        1,1,1,0; 1,1,1,1];
    para.setflagOrder = para.setflag;
    
    para = assignpara(para, varargin{:});
    
    if nnz(para.setflag ~= para.setflagOrder) && para.precalculated
        d1 = numarray2strarray(para.setflag);
        d2 = numarray2strarray(para.setflagOrder);
        k1 = arrayfun(@(x) strjoin(d1(x,:), '_'), 1:size(para.setflag, 1), 'unif', 0);
        k2 = arrayfun(@(x) strjoin(d2(x,:), '_'), 1:size(para.setflag, 1), 'unif', 0);
        [~, si] = ismember(k2, k1);
        para.setflag = para.setflag(si,:);
        setdata = setdata(si);
    end
    
    t = (0:2*pi/39:2*pi);
    angles = [pi/4, pi/3.5, pi*2.5/3.5, pi*3/4];
    offset = [0, 0; -0.8, 0.6; -1.8, 0.6; -2.6, 0];
    labeloffset = [0, 0.2; -0.2, 0.2; -0.4, 0.2; -0.4, 0.2];
    numlabelpos = [0.78,0.46; -0.53,1.63; -0.30,0.74; -2.83,1.63; ...
        -0.67,-1.05; -1.60,0.74; -0.99,-0.16; -3.82,0.46; -1.55,-1.67; ...
        -2.47,-1.05; -1.85,-1.23; -2.90,0.74; -1.15,-1.23; -2.14,-0.16; -1.59,-0.80];    
    
    x = 2.5;
    y = 1;
    
    cla
    hold on;
    %loc = get(gcf, 'position');
    %set(gcf, 'position', [loc(1), loc(2), 750, 650], 'color', 'w');
    handles.patch = NaN(4, 1);
    handles.label = NaN(4, 1);
    handles.numlabel = NaN(15, 1);
    
    for i = 1:length(angles)
        rot = angles(i);
        xd = x*sin(t);
        yd = y*cos(t);
        m = [cos(rot), -sin(rot); sin(rot), cos(rot)];
        xyrot = m*[xd; yd];
        handles.patch(5-i) = patch(xyrot(1,:)+offset(i,1), xyrot(2,:)+offset(i,2), 'k');
        set(handles.patch(5-i), 'facecolor', para.color(i,:), 'facealpha', para.alpha, 'edgecolor', 'none');
        
        handles.label(5-i) = text( xyrot(1,11)+offset(i,1)+labeloffset(i,1), ...
            xyrot(2,11)+offset(i,2)+labeloffset(i,2), para.label{5-i}, 'fontsize', para.labelfontsize);
        twidth = get(handles.label(5-i), 'extent');        
        set(handles.label(5-i), 'position', [xyrot(1,11)+offset(i,1) - twidth(3)/2, ...
            xyrot(2,11)+offset(i,2)+labeloffset(i,2), 0]);        
    end
    
    %fill in numbers of intersections
    setnum.set = para.setflag;
    setnum.num = zeros(size(para.setflag, 1), 1);    
    for i = 1:size(para.setflag,1)
        if ~para.precalculated
            setnum.num(i) = sum(all(bsxfun(@eq, setdata, para.setflag(i,:)), 2));
        else
            setnum.num(i) = setdata(i);
        end
        handles.numlabel(i) = text( numlabelpos(i,1), numlabelpos(i,2), ...
            num2str(setnum.num(i)), 'fontsize', para.numfontsize);
    end
    
    xlim([-4.5 2]);
    ylim([-3 3.5]);
%     set(gca, 'color', 'none')
    axis square 
    axis off
    hold off
end