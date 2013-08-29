function impute = imputeGenotype(data, chrm, pos, postype, theta)
    %Scan through markers and impute missing values by
    %flanking markers and corresponding recombination fraction
    %data: marker x sample; 0 or 1 as genotypes, 2 or NaN as missing
    %chrm: chromosome for each marker
    %pos: position of each marker
    %postype: {'physical','genetic'} for fraction, physical position and 
    %   genetic position (centiMorgan)
    %theta: marker x 1; either recombination fraction (to previous marker)
    %   or physical location; not needed if postype = 'genetic'
    %
    %If theta is recombination fraction, consider conditional
    %probability P(G2=1 | G1, G3, loc1, loc2, loc3).
    %If theta is location(should be genetic location, such as centiMorgan),    
    %resolve G2 by the closer genotype (assuming physical and genetic
    %location are linearly correlated). If G1 and G3 are the same and the
    %distance is small, G2=G1=G3. If G1=G3 but far away, G2 has .5 of
    %probability to be ~G1. But for simplicity, we assume it prefers to be
    %G1.
    %
    
    if nargin < 5, theta = []; end
    
    impute = data;
    [nmarker nsample] = size(impute);
    
    count = 0;
    %indexing markers with missing values and their flanking markers
    missMi = [];
    %#missMi x #sample; each cell is a 2x1 vector [pre_mi, next_mi]        
    flankMi = cell(0,nsample); 
    mask = isnan(impute);
    for mi = 1:nmarker
        if sum(mask(mi,:)) > 0
            missMi = [missMi mi];
            flankMi{end+1,1} = []; %initiate a row
            chrmmask = (chrm == chrm(mi));
            for si = find(mask(mi,:))
                prev_mi = find( ~mask(1:mi-1,si) & chrmmask(1:mi-1), 1, 'last');
                next_mi = find( ~mask(mi+1:end,si) & chrmmask(mi+1:end), 1, 'first');
                if isempty(prev_mi), prev_mi = NaN; end
                if isempty(next_mi)
                    next_mi = NaN; 
                else 
                    next_mi = next_mi + mi;
                end                
                flankMi{end,si} = [prev_mi, next_mi];
            end
        end
    end
        
    if strcmp(postype, 'physical') && ~isempty(theta)
        count = 0; pcount = 0;
        for i = 1:length(missMi)
            for si = 1:nsample
                if isempty(flankMi{i,si}), continue; end  
                li = flankMi{i,si}(1); ri = flankMi{i,si}(2);                              
                p = 1;
                for j = [li ri]
                    if ~isnan(j)
                        frac = calFraction(theta, pos, missMi(i), j);
                        if data(j,si) == 1
                            p = p * (1-frac);
                        else
                            p = p * frac;
                        end
                    end
                end                                
                if sum(isnan(flankMi{i,si})) == 0
                    frac = calFraction(theta, pos, li, ri);                    
                    if data(li,si)==data(ri,si)
                        if 1-frac == 0 
                            if p ~= 0,  count = count + 1;  end                            
                        else
                            p = p / (1-frac);
                        end
                    else
                        if frac == 0 
                            if p ~= 0,  count = count + 1;  end
                        else
                            p = p / frac;
                        end
                    end                    
                end
                if p > 1
                    p = 1;
                    pcount = pcount + 1;
                end
                impute(missMi(i),si) = p;
            end             
        end
        impute(impute>=0.5) = 1;
        impute(impute<0.5) = 0;
        if count > 0
            fprintf('%d values have incorrect recomb fraction info\n',count);
        end
        if pcount > 0
            fprintf('%d values have p>1 \n',pcount);
        end
    elseif strcmp(postype, 'physical')
        for i = 1:length(missMi)
            for si = 1:nsample
                if isempty(flankMi{i,si}), continue; end     
                li = flankMi{i,si}(1); ri = flankMi{i,si}(2);                                           
                if isnan(li)
                    impute(missMi(i),si) = data(ri,si);
                elseif isnan(ri)
                    impute(missMi(i),si) = data(li,si);
                elseif data(li,si) == data(ri,si)
                    impute(missMi(i),si) = data(li,si);
                else
                    if pos(missMi(i))-pos(li) < pos(ri)-pos(missMi(i))
                        impute(missMi(i),si) = data(li,si);
                    else
                        impute(missMi(i),si) = data(ri,si);
                    end
                end
            end 
        end        
    elseif strcmp(postype, 'genetic')
        for i = 1:length(missMi)
            for si = 1:nsample
                if isempty(flankMi{i,si}), continue; end  
                li = flankMi{i,si}(1); ri = flankMi{i,si}(2);                              
                p12 = 1; p23 = 1; p13 = 1;
                if ~isnan(li)
                    p12 = 0.5 + (data(li,si)-0.5) * ...
                        exp(-0.02 * (pos(missMi(i))-pos(li)));
                end
                if ~isnan(ri)
                    p23 = 0.5 + (data(ri,si)-0.5) * ...
                        exp(-0.02 * (pos(ri)-pos(missMi(i))));
                end
                if ~isnan(li+ri)
                    p13 = 0.5 + ((data(li,si)==data(ri,si))-0.5) * ...
                        exp(-0.02 * (pos(ri)-pos(li)));
                end
                impute(missMi(i),si) = p12*p23/p13;
            end 
        end
    else
        error('Unknown theta type %s.\n',thetatype);
    end
end

function frac = calFraction(theta, pos, si, ei)
    if ei < si
        tmp = si;
        si = ei;
        ei = tmp;
    end
    frac = 0;
    for i = si+1:ei
        frac = frac + theta(i)*(pos(i)-pos(i-1)+1);
    end
    frac = frac / (pos(ei)-pos(si)+1);
end
    