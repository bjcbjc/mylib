function res = extractpheno(data, time, varargin)
    %data: #timept x # series
    
    para.swindow = 15;
    para.include = 0.20; %the max slope +- include*ratio
    para.relaxinclude = 0.30;
    para.expand = 5; %expand 5 pts each iteration (one direction)
    para.expanditer = 10; %max #iteration to expand in one direction
    para.minnpt = 20;
    para.exclude = [10 50]; %exclude the first # and the end # data pts
    para.blank = 0;
    para.validrange = [0.05 0.5]; % 6~7 doubling (after substraction from the min or blank)
    para.adaptrange = 0.15;
    para.adaptmaxgap = 3;
    
    para = assignpara(para, varargin{:});
    
    %remove NaN chunk from data first
    i = find(isnan(time), 1, 'first');
    time(i:end) = [];
    data(i:end,:) = [];
    
    [npt ndata] = size(data);
    assert(size(time,1)==npt, '#time pt should be the same\n');
    
    para.QCi_relaxwindow = 1;
    para.QCi_chooseLong = 2;
    para.QCi_shortSeg = 3;
    para.QCi_invalid = 4;
    
    [pdata datamask] = preprocess(data, para);
    
    res.maxslope = NaN(ndata,1);
    res.peakpti = NaN(ndata,1);
    res.pti = NaN(ndata, 2);
    res.regb = NaN(ndata, 2);
    res.QC = false(ndata, 4); %
    res.longpti = NaN(ndata, 2);
    res.sw = NaN(ndata, 1);
    res.adapt = NaN(ndata, 1);
    res.para = para;
    
    for di = 1:ndata
        %get valid data
        vi = find(datamask(:, di) == 1);
        
        if length(vi) < para.minnpt
            res.QC(di, para.QCi_invalid) = true;
            continue
        end
        
                
        cursw = para.swindow;
        
        while 1
            sdata = log2(smooth(pdata(:, di), cursw));
            allslope = smooth(gradient(sdata, time), 5);
            slope = allslope(vi); %because we include two end points in vdata
            
            [maxslope peakpti] = max(slope);
            %res.peakpti(di) = res.peakpti(di) + para.exclude(1);
            diffratio = (slope - maxslope) / maxslope;

            [peaksegpti longpti longseg] = expsegment(diffratio, para.include, peakpti);
            if peaksegpti(2)-peaksegpti(1)+1 < para.minnpt
                res.QC(di, para.QCi_relaxwindow) = true;
                [peaksegpti longpti longseg] = expsegment(diffratio, para.relaxinclude, peakpti);
                if peaksegpti(2) - peaksegpti(1) + 1 < para.minnpt
                    if cursw >= 23
                        break
                    end
                    cursw = cursw + 2;
                else
                    break
                end
%                 %still very short, force to expand at least half of the
%                 %minnpt
%                 curnpt = peaksegpti(2)-peaksegpti(1)+1;
%                 peaksegpti(1) = max(1, peaksegpti(1) - ceil((para.minnpt-curnpt)/2));
%                 peaksegpti(2) = min(npt, peaksegpti(2) + ceil((para.minnpt-curnpt)/2));
%             end
            else
                break
            end
        end
        
        %the longest segment and the peak segment are different
        %look to see which segment gives average higher slope (in case the
        %data is noisy)       
        %adjust the window and fit the curve
        if longseg
            [finalpti regb] = iterfitting(...
                sdata(vi), time(vi), peaksegpti(1):peaksegpti(2), para, peakpti);
        else 
            avgpeak = mean(slope(setdiff(peaksegpti(1):peaksegpti(2), peakpti)));
            tmp = sort(slope(longpti(1):longpti(2)));
            avglong = mean(tmp(2:end));
            if avgpeak >= avglong
                [finalpti regb] = iterfitting(...
                    sdata(vi), time(vi), peaksegpti(1):peaksegpti(2), para, peakpti);
            else
                [finalpti regb] = iterfitting(...
                    sdata(vi), time(vi), longpti(1):longpti(2), para, peakpti);
                res.QC(di, para.QCi_chooseLong) = true;
            end                
        end
            
        
        %put the result into res
        ptishift = vi(1) - 1;
        res.maxslope(di) = maxslope;
        res.peakpti(di) = peakpti + ptishift;
        res.pti(di,:) = finalpti + ptishift;
        res.regb(di,:) = regb;
        res.longpti(di,:) = longpti + ptishift;
        res.sw(di) = cursw;
        res.adapt(di) = adapttime(allslope(1:res.pti(di,1)), regb(2), para);
        if res.pti(di,2)-res.pti(di,1)+1 < para.minnpt
            res.QC(di, para.QCi_shortSeg) = true;
        end
    end
    
end

function [pdata datamask] = preprocess(data, para)    
    [n d] = size(data);
    pdata = data - para.blank; 
    pdata = pdata - repmat(min(pdata(1:n-para.exclude(2),:)), n, 1); %avoid the weird "-1" in the end
    %substraction of the minimum is used to constrain the data used is
    %within 6~7 doublings (0.05~0.5); but calculation of DT needs the
    %original OD
        
    datamask = true(n, d);
    datamask(1:para.exclude(1), :) = false;
    %datamask(n-para.exclude(2)-1, :) = false; %before 100303
    datamask(n-para.exclude(2)-1:end, :) = false;
    datamask = datamask & (pdata >= para.validrange(1)) & (pdata <= para.validrange(2));
    
    for i = 1:d
        datamask(find(data(:,i)>0.7, 1, 'first'):end, i) = false;
    end
    
    %correct to have continuous valid range
    for i = 1:d
        zerotick = find(~datamask(:,i));
        [tmp mi] = max(diff(zerotick));
        si = zerotick(mi)+1;
        ei = zerotick(mi+1)-1;
        assert(ei-si+1==tmp-1, 'error finding valid range');
        datamask(1:si-1, i) = false;
        datamask(ei+1:end, i) = false;
    end    
    
    pdata = data - min( repmat(min(data(1:n-para.exclude(2),:)), n, 1), para.blank); 
    %we don't need to shift the min of each curve to zero (as in
    %pdata=pdad-min(pdata), because that will affect the calculation of DT;
    %but here we don't want negative OD, do check if the curve has minimum
    %value smaller than blank; if there is, substract the data from its
    %minimum on

end

function adpt = adapttime(slope, growthrate, para)
    ratio = (slope - growthrate) / growthrate;
    pti = find(ratio >= -para.adaptrange & ratio <= para.adaptrange);
    if length(pti) < 2
        adpt = length(slope); %take the start of the exp phase
    else
        endofgap = find(diff(pti) > para.adaptmaxgap, 1, 'last');
        if ~isempty(endofgap)
            adpt = pti(endofgap+1);
        else
            adpt = pti(1);
        end
    end
end

function [peaksegpti longpti longseg] = expsegment(diffratio, ratiocut, peakpti)
    %potential points to be included
    longseg = true;
    pti = find(diffratio >= -ratiocut & diffratio <= ratiocut);
    peaksegpti = [pti(1) pti(end)];
    longpti = [pti(1) pti(end)];
    %if the segment is not continuous
    if any(diff(pti)>1)
        segstartindex = [1; find( diff(pti) > 1)+1 ];
        setmaxpti = find( pti == peakpti );
        %determine the length of each segment (each segment is defined as a
        %continuous region that are within the defined ratio of max slope)
        seglen = [diff( segstartindex ); length(pti)-segstartindex(end)];

        %find the longest segment
        [maxseglen segi] = max(seglen);        
        if segi == length(seglen) %the last segment
            longpti = [pti(segstartindex(end)) pti(end)];       
        else %the longest segment is the first or in the middle
            longpti = [pti(segstartindex(segi)) pti(segstartindex(segi+1)-1)];
        end

        %include the peakpti in the data
        peaksegi = find(segstartindex <= setmaxpti, 1, 'last');

        if peaksegi == length(segstartindex)
            peaksegpti = [pti(segstartindex(end)) pti(end)]; %last segment
        else
            peaksegpti = [pti(segstartindex(peaksegi)) pti(segstartindex(peaksegi+1)-1)];
        end
        if max([longpti(1) peaksegpti(1)]) > min([longpti(2) peaksegpti(2)])
            longseg = false; %the selected segment is not the longest segment
        end
    end
    %peaksegpti(2) = peaksegpti(2) + 1;
    %longpti(2) = longpti(2) + 1;
end

function [bestpti, bestregb, minerr] = iterfitting(sdata, time, pti, para, peakpti)
    %iteratively expand/shrink the points to fit the exponential curve
    minerr = Inf;
    npt = length(sdata);
    for iterpti = -para.expanditer*para.expand:para.expand:para.expanditer*para.expand
        for iterptj = -para.expanditer*para.expand:para.expand:para.expanditer*para.expand
            si = max(1, pti(1)+iterpti);
            ei = min(npt, pti(end)+iterptj);
            if si >= ei || (si>pti(1) || ei<pti(end) && ei-si<para.minnpt) ...
                    || (ismember(peakpti, pti) && (si>peakpti || ei<peakpti))
                continue
            end
            [regb tmp tmp tmp st] = regress(sdata(si:ei), [ones(ei-si+1,1) time(si:ei)]);
%            residual = 2.^sdata(si:ei) - 2.^([ones(ei-si+1,1) time(si:ei)]*regb);
%             r2 = corr(sdata(si:ei), [ones(ei-si+1,1) time(si:ei)]*regb)^2;
            %residual = sdata(si:ei) - [ones(ei-si+1,1) time(si:ei)]*regb;
            %like = neglikecost(residual);
%            if mean(residual.^2) < minerr
%              if -r2 < minerr
             if st(3) < minerr
%                minerr = mean(residual.^2);
%                 minerr = -r2;
                minerr = st(3);
                bestpti = [si ei];
                bestregb = regb;
            end
        end
    end
end