function fhandles = displayData3(genotype, exp, growth, fhandles, climexp, climgrowth)
    %plot data on three windows
    %fhandles: [fgeno, fexp, fgrowth]; [] will generate new fhandles
    %climexp: clim for exp data
    %climgrowth: clim for growth
    %
    
    if nargin < 6
        climgrowth = [];
    end
    if nargin < 5
        climexp = [];
    end
    if nargin < 4
        fhandles = [];
    end
    if isempty(fhandles)
        for i = 1:3
            fhandles(i) = figure();
        end
    end
    figure(fhandles(1))
    nanimagesc(genotype, genColorMap('bw',2));
    
    figure(fhandles(2))
    if isempty(climexp)
        nanimagesc(exp, genColorMap('rkg',128));        
    else
        nanimagesc(exp, genColorMap('rkg',128),'clim',climexp);        
    end
    colorbar;
    
    figure(fhandles(3))
    if isempty(climgrowth)
        nanimagesc(growth, genColorMap('yb',128));
    else
        nanimagesc(growth, genColorMap('yb',128));
    end
    colorbar;
end