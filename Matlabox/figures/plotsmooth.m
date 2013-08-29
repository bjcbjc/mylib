function plotsmooth(varargin)

    para.pt = 1000;
    para.method = 'cubic';
    para.passpara = {};
    
    if nargin == 1
        y = varargin{1};
        x = (1:length(y))';
    else
        x = varargin{1};
        y = varargin{2};
        for i = 3:2:length(varargin)
            if isfield(para, varargin{i})
                para.(varargin{i}) = varargin{i+1};            
            else
                if i+1 <= length(varargin)
                    para.passpara = [para.passpara; varargin{i}; varargin{i+1}];
                else
                    para.passpara = [para.passpara; varargin{i}];
                end
            end
        end
    end
        
    [nx ny] = size(x);
    [nx2 ny2] = size(y);
    
    if isvector(x) && nx == 1
        x = x';
        nx = ny;
        ny = 1;
    end
    if isvector(y) && nx2 == 1
        y = y';
        nx2 = ny2;
        ny2 = 1;
    end
    
    assert(nx==nx2, 'x and y must have the same number of samples');
    if ny > 1
        assert(ny==ny2, 'x should have one column or the same number of column as y');
    end
        
    if ny == 1
        plotx = (x(1) : (x(end)-x(1))/para.pt : x(end))';
        ploty = interp1(x, y, plotx, 'cubic');
    else
        for yi = 1:ny
            plotx(:,yi) = (x(1,yi) : (x(end,yi)-x(1,yi))/para.pt : x(end,yi))';
            ploty(:,yi) = interp1(x(:,yi), y(:,yi), plotx(:,yi), 'cubic');
        end
    end
    plot(plotx, ploty, para.passpara{:});

end