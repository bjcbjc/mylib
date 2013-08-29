function [chrms chrmtick] = plotgenomescore(score, pos, varargin)
    %score: #test x #score
    %pos: #test x 2, positions of markers used for test; first column is
    %   chrm, and the second is position in the chrm
    %varargin: 
    %   .peak: # of peak to note the position, def = 1;
    %   .peakdist: min dist between peaks, def = 200k;
    %   .xlablel, .ylabel
    %

    [nm ns] = size(score);
    assert(nm==size(pos,1), '#maker does not match the dimension of score');
    
    para.peak = 1;
    para.peakdist = 200000;
    para.xlabel = 'chrm';
    para.ylabel = '-log10(p)';
    para.plotfun = @plot;
    para.xlim = [];
    para.ylim = [];
    
    para = assignpara(para, varargin{:});

    [pos rowi] = sortrows(pos);
    score = score(rowi, :);
    chrms = unique(pos(:,1));
    nchrm = length(chrms);
    
    chrmtick = NaN(nchrm,2);
    for i = 1:nchrm
        chrmtick(i,1) = find(pos(:,1)==chrms(i), 1, 'first');
        chrmtick(i,2) = find(pos(:,1)==chrms(i), 1, 'last');
    end
        
    for i = 1:ns
        subplot(ns,1,i);
        para.plotfun(score(:,i), 'x');
        if ~isempty(para.xlim)
            xlim(para.xlim);
        end
        if ~isempty(para.ylim)
            ylim(para.ylim);
        end
        ym = ylim();
        xm = xlim();
        M = max(nanmax(score(:,i)),ym(2))+1;
        hold on
        for li = 1:nchrm-1
            plot([chrmtick(li,2)-0.5 chrmtick(li,2)-0.5],[0 M],'--','color',[0.7 0.7 0.7]);
        end
        pindexes = findpeak(score(:,i), para.peak, para.peakdist, pos);
        for pi = 1:length(pindexes)
            if ~isnan(pindexes(pi))
                para.plotfun(pindexes(pi),score(pindexes(pi),i),'rx');
                text(pindexes(pi)+(xm(2)-xm(1))/50, ...
                    score(pindexes(pi),i)+(ym(2)-ym(1))/50, ...
                    sprintf('%.1fK',pos(pindexes(pi),2)/1000));
            end
        end
        hold off
        set(gca, 'xtick', mean(chrmtick,2), 'xticklabel', chrms);
        xlabel(para.xlabel);
        ylabel(para.ylabel);
    end
end


function indexes = findpeak(s, np, pdist, pos)
    c = 0;
    indexes = NaN(np,1);
    while c < np && c < length(s)
        c = c + 1;
        [tmp indexes(c)] = nanmax(s);
        s(pos(:,1)==pos(indexes(c),1) & ...
            pos(:,2) >= pos(indexes(c),2) - pdist/2 & ...
            pos(:,2) <= pos(indexes(c),2) + pdist/2) = NaN;
    end
end