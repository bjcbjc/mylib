function status = savePdfFromEps(fig, fn, keepEps)
    if nargin < 3
        keepEps = true;
    end
    if strcmp(fn(end-3:end), '.pdf')
        epsFn = [fn(1:end-4), '.eps'];
    elseif strcmp(fn(end-3:end), '.eps')
        epsFn = fn;
        fn = [fn(1:end-4), '.pdf'];
    else
        epsFn = [fn, '.eps'];
        fn = [fn, '.pdf'];
    end
                
    saveas(fig, epsFn, 'epsc2');
    cmd = 'gs -dEPSCrop -dBATCH -dSAFER -dNOPAUSE ';    
    cmd = sprintf('%s -dAutoRotatePages=/None -q -sDEVICE=pdfwrite -sOutputFile=%s %s', cmd, fn, epsFn);
    status = system(cmd);
    if ~keepEps
        system(sprintf('rm -f %s', epsFn));
    end
end