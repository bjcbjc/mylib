function [handles, setnum] = fixvenn(setdata, varargin)
    %logicaldata: #data x #set, should be logical (or will be converted to
    %   logical)
    %
    %
    
    if ~islogical(setdata)
        setdata = setdata ~= 0;
    end
    nset = size(setdata,2);
    
    
    %para = assignpara(para, varargin{:});
    switch nset
        case 4
            [handles, setnum] = venn4(setdata, varargin{:});
        otherwise
            fprintf('%d sets are not supported\n', nset);
    end
    
end


function [handles, setnum] = venn4(setdata, varargin)
    para.label = repmat({''}, 1, 4);
    para.color = [0.7, 0.3, 0.4; 0.3, 0.7, 0.4; 0.7, 0.6, 0.3; 0.4, 0.3, 0.7];
    para.alpha = 0.5;
    para.numfontsize = 16;
    para.labelfontsize = 16;
    
    para = assignpara(para, varargin{:});
    
    t = (0:2*pi/39:2*pi);
    angles = [pi/4, pi/3.5, pi*2.5/3.5, pi*3/4];
    offset = [0, 0; -0.8, 0.6; -1.8, 0.6; -2.6, 0];
    labeloffset = [0, 0.2; -0.2, 0.2; -0.4, 0.2; -0.4, 0.2];
    numlabelpos = [0.78,0.46; -0.53,1.63; -0.30,0.74; -2.83,1.63; ...
        -0.67,-1.05; -1.60,0.74; -0.99,-0.16; -3.82,0.46; -1.55,-1.67; ...
        -2.47,-1.05; -1.85,-1.23; -2.90,0.74; -1.15,-1.23; -2.14,-0.16; -1.59,-0.80];
    setflag = [0,0,0,1; 0,0,1,0; 0,0,1,1; 0,1,0,0; 0,1,0,1; 0,1,1,0; ...
        0,1,1,1; 1,0,0,0; 1,0,0,1; 1,0,1,0; 1,0,1,1; 1,1,0,0; 1,1,0,1; ...
        1,1,1,0; 1,1,1,1];
    
    x = 2.5;
    y = 1;
    
    clf
    hold on;
    loc = get(gcf, 'position');
    set(gcf, 'position', [loc(1), loc(2), 750, 650]);
    handles.patch = NaN(4, 1);
    handles.label = NaN(4, 1);
    handles.numlabel = NaN(15, 1);
    
    for i = 1:length(angles)
        rot = angles(i);
        xd = x*sin(t);
        yd = y*cos(t);
        m = [cos(rot), -sin(rot); sin(rot), cos(rot)];
        xyrot = m*[xd; yd];
        handles.patch(i) = patch(xyrot(1,:)+offset(i,1), xyrot(2,:)+offset(i,2), 'k');
        set(handles.patch(i), 'facecolor', para.color(i,:), 'facealpha', para.alpha, 'edgecolor', 'none');
        
        handles.label{i} = text( xyrot(1,11)+offset(i,1)+labeloffset(i,1), ...
            xyrot(2,11)+offset(i,2)+labeloffset(i,2), para.label{i}, 'fontsize', para.labelfontsize);
        twidth = get(handles.label{i}, 'extent');        
        set(handles.label{i}, 'position', [xyrot(1,11)+offset(i,1) - twidth(3)/2, ...
            xyrot(2,11)+offset(i,2)+labeloffset(i,2), 0]);        
    end
    
    %fill in numbers of intersections
    setnum.set = setflag;
    setnum.num = zeros(size(setflag, 1), 1);    
    for i = 1:size(setflag,1)
        setnum.num(i) = sum(all(bsxfun(@eq, setdata, setflag(i,:)), 2));
        handles.numlabel{i} = text( numlabelpos(i,1), numlabelpos(i,2), ...
            num2str(setnum.num(i)), 'fontsize', para.numfontsize);
    end
    
    xlim([-4.5 2.5]);
    ylim([-2.5 4.5]);
    set(gca, 'color', 'none')
    axis square 
    axis off
    hold off
end