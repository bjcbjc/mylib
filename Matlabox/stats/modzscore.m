function modz = modzscore(obsData, sampleData, dim)
    %calcualte modified z-score of obsData to sampleData (distribution)
    %
    % if sampleData is a matrix, dim suggests which dimension is "variable"
    % and which is samples of the variable
    % if obsData is a vector or a matrix, and sampleData is a matrix, the
    % "variable" dimension of obsData much match sampleData:
    %   obsData: variable x obs, sampleData: variable x samples, modz:
    %       variable x obs
    %   obsData: 1 x obs, sampleData: variable x sample, modz: 1 x obs
    %   
    %
    
    medAbsDev = mad(sampleData, 1, dim) * 1.486; %median abs dev
    meanAbsDev = mad(sampleData, 0, dim) * 1.253314;
    i = medAbsDev == 0;
    medAbsDev(i) = meanAbsDev(i);
    
    med = nanmedian(sampleData, dim);
    
    tmp = bsxfun(@minus, obsData, med);
    modz = bsxfun(@rdivide, tmp, medAbsDev);
    
    if any(medAbsDev == 0)
        infidx = bsxfun(@and, medAbsDev==0, tmp > 0);
        modz(infidx) = Inf;
        infidx = bsxfun(@and, medAbsDev==0, tmp < 0);
        modz(infidx) = -Inf;
    end
    