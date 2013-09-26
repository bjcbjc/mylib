function writeHeatmapImage(mtx, fn, varargin )
    para.width = 20;
    para.height = 12;
    para.fontsize = 10;
    para.title = '';
    para.rowlabel = {};
    para.collabel = {};
    para.colormap = redbluecmap(20);
    para.showcolorbar = true;
    para.clim = [];
    para.xmargin = [120, 50]; %left, right
    para.ymargin = [100, 50]; %top, bottom
    para.grid = true;

    para = assignpara(para, varargin{:});
    
    
    if isempty(para.clim)
        para.clim = prctile(mtx(~isnan(mtx)), [5, 95]);
    end
        
    h = figure('visible', 'off');
    [m, n] = size(mtx);
    set(h, 'position', [0, 0, n*para.width+sum(para.xmargin), m*para.height+sum(para.ymargin)], 'paperpositionmode', 'auto');
    set(gca, 'units', 'pixels', 'position', [para.xmargin(1), para.ymargin(2), n*para.width, m*para.height], ...
        'color', [0.8, 0.8, 0.8]);
    
    
    nanimagesc(mtx, para.colormap, 'clim', para.clim);
    if para.grid
        gridforimage('color', [0.7,0.7,0.7]);
    end
    if para.showcolorbar
%         cbarh = colorbar( 'units', 'pixels', 'position', ...
%             [para.xmargin(1)+(n+1)*para.width, ...
%             m*para.height+sum(para.ymargin) - para.ymargin(1)-para.height*10, ...
%             1, ...
%             para.height*10]);
        cbarh = colorbar( 'location', 'east');
        set(cbarh, 'unit', 'pixels', 'position', ...
            [para.xmargin(1)+(n+1)*para.width, ...
            m*para.height+sum(para.ymargin) - para.ymargin(1)-para.height*10, ...
            10, ...
            para.height*10]);
    end
    
    if ~isempty(para.rowlabel)
        set(gca, 'ytick', 1:m, 'yticklabel', para.rowlabel, 'fontsize', para.fontsize);
    end
    if ~isempty(para.collabel)
        texth = topVertXlabel(para.collabel, 90, 0);
    else
        texth = [];
    end
    if ~isempty(para.title)
        if isempty(texth)
            title(para.title, 'fontsize', para.fontsize, 'fontweight', 'bold');
        else
            ypos = cell2mat(get(texth,'extent'));
            ypos = -max(ypos(:,4))-1;
            text(n/2, ypos, para.title, 'fontsize', para.fontsize, 'fontweight', 'bold','unit','data');
        end
    end
    saveas(h, fn, 'png');
    close(h);    
end

function hText = topVertXlabel(xTickLabels,rot,y, varargin)
    %hText = xticklabel_rotate(XTick,rot,XTickLabel,varargin)     Rotate XTickLabel
    %
    % Syntax: xticklabel_rotate
    %
    % Input:    
    % {opt}     XTick       - vector array of XTick positions & values (numeric) 
    %                           uses current XTick values or XTickLabel cell array by
    %                           default (if empty) 
    % {opt}     rot         - angle of rotation in degrees, 90° by default
    % {opt}     XTickLabel  - cell array of label strings
    % {opt}     [var]       - "Property-value" pairs passed to text generator
    %                           ex: 'interpreter','none'
    %                               'Color','m','Fontweight','bold'
    %
    % Output:   hText       - handle vector to text labels
    %    
    %
    % Note : you can not RE-RUN xticklabel_rotate on the same graph. 
    %
    % This is a modified version of xticklabel_rotate90 by Denis Gilbert
    % Modifications include Text labels (in the form of cell array)
    %                       Arbitrary angle rotation
    %                       Output of text handles
    %                       Resizing of axes and title/xlabel/ylabel positions to maintain same overall size 
    %                           and keep text on plot
    %                           (handles small window resizing after, but not well due to proportional placement with 
    %                           fixed font size. To fix this would require a serious resize function)
    %                       Uses current XTick by default
    %                       Uses current XTickLabel is different from XTick values (meaning has been already defined)
    % Brian FG Katz
    % bfgkatz@hotmail.com
    % 23-05-03
    % Modified 03-11-06 after user comment
    %	Allow for exisiting XTickLabel cell array
    % Modified 03-03-2006 
    %   Allow for labels top located (after user comment)
    %   Allow case for single XTickLabelName (after user comment)
    %   Reduced the degree of resizing
    % Modified 11-jun-2010
    %   Response to numerous suggestions on MatlabCentral to improve certain
    %   errors.
    % Other m-files required: cell2mat
    % Subfunctions: none
    % MAT-files required: none
    %
    % See also: xticklabel_rotate90, TEXT,  SET
    % Based on xticklabel_rotate90
    %   Author: Denis Gilbert, Ph.D., physical oceanography
    %   Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
    %   email: gilbertd@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
    %   February 1998; Last revision: 24-Mar-2003
    %
    % Modified by BJ Chen 2013
    %
   
    if nargin < 2
        rot = 90 ;
    end

    %Make XTick a column vector
    XTick = (1:length(xTickLabels))';

    %Set the Xtick locations and set XTicklabel to an empty string
    set(gca,'XTick',XTick,'XTickLabel','')

    
    % Determine the location of the labels based on the position
    % of the xlabel
    hxLabel = get(gca,'XLabel');  % Handle to xlabel
    set(hxLabel,'Units','data');
%     xLabelPosition = get(hxLabel,'Position');
%     y = xLabelPosition(2);

    %CODE below was modified following suggestions from Urs Schwarz
    y=repmat(y,size(XTick,1),1);
    % retrieve current axis' fontsize
    fs = get(gca,'fontsize');

    % Place the new xTickLabels by creating TEXT objects
    hText = text(XTick, y, xTickLabels,'fontsize',fs);

    % Rotate the text objects by ROT degrees
    %set(hText,'Rotation',rot,'HorizontalAlignment','right',varargin{:})
    % Modified with modified forum comment by "Korey Y" to deal with labels at top
    % Further edits added for axis position
%     xAxisLocation = get(gca, 'XAxisLocation');  
% 
%     if strcmp(xAxisLocation,'bottom')  
%         set(hText,'Rotation',rot,'HorizontalAlignment','right',varargin{:})  
%     else  
        set(hText,'Rotation',rot,'HorizontalAlignment','left',varargin{:})  
%     end

    % Adjust the size of the axis to accomodate for longest label (like if they are text ones)
    % This approach keeps the top of the graph at the same place and tries to keep xlabel at the same place
    % This approach keeps the right side of the graph at the same place 

    set(get(gca,'xlabel'),'units','data')           ;
        labxorigpos_data = get(get(gca,'xlabel'),'position')  ;
    set(get(gca,'ylabel'),'units','data')           ;
        labyorigpos_data = get(get(gca,'ylabel'),'position')  ;
    set(get(gca,'title'),'units','data')           ;
        labtorigpos_data = get(get(gca,'title'),'position')  ;

    %set(gca,'units','pixel')                        ;
    set(hText,'units','pixel')                      ;
    set(get(gca,'xlabel'),'units','pixel')          ;
    set(get(gca,'ylabel'),'units','pixel')          ;

    origpos = get(gca,'position')                   ;

    % textsizes = cell2mat(get(hText,'extent'))       ;
    % Modified with forum comment from "Peter Pan" to deal with case when only one XTickLabelName is given. 
    x = get( hText, 'extent' );  
    if iscell( x ) == true  
        textsizes = cell2mat( x ) ;  
    else  
        textsizes = x;  
    end  

%     largest =  max(textsizes(:,3))                  ;
    longest =  max(textsizes(:,4))                  ;
    
%     laborigext = get(get(gca,'xlabel'),'extent')    ;
%     laborigpos = get(get(gca,'xlabel'),'position')  ;
% 
%     labyorigext = get(get(gca,'ylabel'),'extent')   ;
%     labyorigpos = get(get(gca,'ylabel'),'position') ;
%     leftlabdist = labyorigpos(1) + labyorigext(1)   ;

    % assume first entry is the farthest left
    leftpos = get(hText(1),'position')              ;
    leftext = get(hText(1),'extent')                ;
    leftdist = leftpos(1) + leftext(1)              ;
    if leftdist > 0,    leftdist = 0 ; end          % only correct for off screen problems

    % botdist = origpos(2) + laborigpos(2)            ;
    % newpos = [origpos(1)-leftdist longest+botdist origpos(3)+leftdist origpos(4)-longest+origpos(2)-botdist]  
    %
    % Modified to allow for top axis labels and to minimize axis resizing
%     if strcmp(xAxisLocation,'bottom')     
%         newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
%             origpos(2)+((longest+laborigpos(2))-get(gca,'FontSize')) ...
%             origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
%             origpos(4)-((longest+laborigpos(2))-get(gca,'FontSize'))] ;    
%     else
%         newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
%                 origpos(2) ...
%                 origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
%                 origpos(4)-(longest)+get(gca,'FontSize')]  ;
%     end
    %set(gca,'position',newpos)                      ;

    % readjust position of text labels after resize of plot
    set(hText,'units','data')                       ;
    for loop= 1:length(hText),
        set(hText(loop),'position',[XTick(loop), y(loop)])  ;
    end

    % adjust position of xlabel and ylabel
    laborigpos = get(get(gca,'xlabel'),'position')  ;
    set(get(gca,'xlabel'),'position',[laborigpos(1) laborigpos(2)-longest 0])   ;

    % switch to data coord and fix it all
    set(get(gca,'ylabel'),'units','data')                   ;
    set(get(gca,'ylabel'),'position',labyorigpos_data)      ;
    set(get(gca,'title'),'position',labtorigpos_data)       ;

    set(get(gca,'xlabel'),'units','data')                   ;
        labxorigpos_data_new = get(get(gca,'xlabel'),'position')  ;
    set(get(gca,'xlabel'),'position',[labxorigpos_data(1) labxorigpos_data_new(2)])   ;


    if nargout < 1,
        clear hText
    end
end