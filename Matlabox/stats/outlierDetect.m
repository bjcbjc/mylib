function [out_index] = outlierDetect(vec, varargin)

% NOTE: vec should NOT have NaN
% revised: vec can contain NaN

    valid_index = find(isnan(vec)==0);
    out_index = valid_index;
    i = 1;
    while i <= length(varargin)
        switch lower(varargin{i})
            case 'sigma'
                [nrm, sidx] = sort(zeromean_univar_normalization(vec(valid_index)));
                out = valid_index(sidx(sigma3Outlier(nrm)));
            case 'nrm1'
                [nrm, sidx] = sort(zeromean_univar_normalization(vec(valid_index)));
                out = valid_index(sidx(diffOutlier1(nrm, 1, 2)));
            case 'srt2'
                [srt, sidx] = sort(vec(valid_index));
                out = valid_index(sidx(diffOutlier2(srt, 2)));
            case 'srt1'
                [srt, sidx] = sort(vec(valid_index));
                %out = valid_index(sidx(diffOutlier3(srt, 0.09, 2)));
                %out = valid_index(sidx(diffOutlier3(srt, 0.2, 2)));                
                out = valid_index(sidx(diffOutlier3(srt, std(srt)*2, 2)));
            case 'iqr'
                out = valid_index(iqrOutlier(vec(valid_index)));
        end
        out_index = intersect(out_index, out);
        i = i + 1;
    end        
end

function [outidx] = sigma3Outlier(vec) 
%vec: sorted, normalized vector (zero_mean, uni_var)
%filter out values > 3*sigma
    outidx = find(abs(vec)>3);
end

function [outidx] = diffOutlier1(vec, thres, ntail)
%vec: sorted, normalized vector
%thres: threshold for 1st-order-diff cutoff
%ntail: # of points consider at both ends
%detect abrupt increase (decrease)
    dif = diff(vec);   
    outidx = [];
    for i = ntail:-1:1        
        if dif(i) >= 1.5*thres 
            outidx = [outidx 1:i];
            break
        elseif dif(i) >= thres
            if i == 1 || dif(i-1) < thres
                outidx = [outidx 1:i];
                break
            end
        end
    end    
    len = length(vec);
    for i = len-ntail:1:len-1
        if dif(i) >= 1.5*thres 
            outidx = [outidx i+1:len];
            break
        elseif dif(i) >= thres
            if i == len-1 || dif(i+1) < thres
                outidx = [outidx i+1:len];
                break
            end
        end
    end    
end

function [outidx] = diffOutlier2(vec, ntail)
%vec: sorted vector
%thres: threshold for 1st-order-diff cutoff
%ntail: # of points consider at both ends
%detect abrupt increase (decrease) in the 2nd order diff
    dif = abs(diff(vec));   
    dif2 = abs(diff(dif));
    outidx = [];
    for i = ntail:-1:1        
        if dif2(i) >= 2*dif2(i+1)
            outidx = [outidx 1:i];
            break        
        end
    end    
    len = length(vec);
    for i = len-ntail:1:len-1
        if dif2(i-1) >= 2*dif2(i-2)
            outidx = [outidx i+1:len];
            break        
        end
    end    
end

function [outidx] = diffOutlier3(vec, thres, ntail)
%vec: sorted vector
%thres: threshold for 1st-order-diff cutoff
%ntail: # of points consider at both ends
%detect abrupt increase (decrease)
    dif = abs(diff(vec));       
    outidx = [];
    for i = ntail:-1:1             
        if dif(i) > thres            
            outidx = [outidx 1:i];
            break        
        end
    end    
    len = length(vec);
    for i = len-ntail:1:len-1        
        if dif(i) > thres
            outidx = [outidx i+1:len];
            break        
        end
    end    
end

function [outidx] = iqrOutlier(vec)
%iqr, prctile ignore NaN directly

    qr = 3*iqr(vec);
    Q1 = prctile(vec,25);
    Q3 = prctile(vec,75);
    outidx = [];
    for i = 1:length(vec)
        if vec(i) < (Q1-qr) || vec(i) > (Q3+qr)
            outidx = [outidx i];
        end
    end
end