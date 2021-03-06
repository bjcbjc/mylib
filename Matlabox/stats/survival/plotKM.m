function [curvehandle, pval] = plotKM(survivalData, sample, alteredSample, nonAlteredSample, label, idx)
    if nargin < 6
        idx = false;
    end
    
    if ~idx
        validSample = ~isnan(survivalData(:,1));
        alterIdx = ismember(sample,alteredSample) & validSample;
        nonAlterIdx = ismember(sample,nonAlteredSample) & validSample;
    else
        alterIdx = alteredSample;
        nonAlterIdx = nonAlteredSample;
    end
    [~, ~, ~, ~, plotdata] = logrank( ...
        survivalData(alterIdx,:), ...
        survivalData(nonAlterIdx,:), ...
        0.05, 0, 0);    
    
    curvehandle = flexKmPlot(plotdata.time, plotdata.surv, plotdata.censorXY, ...
        'linecolor',{'r','b'}, 'censorcolor',{'r','b'},...
        'censormark', '+', ...
        'legend', ...
        {sprintf('altered (n=%d)',sum(alterIdx)), ...
        sprintf('non-altered (n=%d)', sum(nonAlterIdx)) } );
    pval = mylogrank(survivalData(alterIdx,:), survivalData(nonAlterIdx,:));
    xlabel('month', 'fontsize',12);
    ylabel(label, 'fontsize', 12);
    title(sprintf('logrank pval=%g',pval),'fontsize',12);