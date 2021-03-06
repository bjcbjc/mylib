classdef batchEffectPCA < handle

    methods (Static)
        function report = test(readCount, covariate, varargin)
            % readCount: #gene x #sample
            % covariate: #sample x #variable            
            para.standardize = true;
            para.log2transform = true;
            para.test = @anova1;
            para.nPC = 20;
            para.covLabel = {};
            para.returnWeight = false;
                        
            para = assignpara(para, varargin{:});            
            para.nPC = min( para.nPC, size(readCount,2)-1);
            
            if ischar(para.test)
                para.test = str2func(para.test);
            end
            
            if para.log2transform
                readCount = log2( readCount + 1);
            end
            if para.standardize
                readCount = zeromean_univar_normalization(readCount, 2);
            end
            report.validGeneIdx = ~any(isnan(readCount),2);
            if para.returnWeight
                [report.pc, sc, l] = pca(readCount(report.validGeneIdx,:)');
            else
                [~, sc, l] = pca(readCount(report.validGeneIdx,:)');
            end
%             [~, sc, l] = princomp(readCount(report.validGeneIdx,:)', 'econ');            
            l = l ./ sum(l);
            report.para = para;
            report.covariate = covariate;
            report.pcVarExplained = l(1:para.nPC);
            report.pcProjection = sc(:, 1:para.nPC);            
            nvar = size(covariate, 2);
            
            if strcmp(func2str(para.test), 'corr')
                [~, report.pcCovAssociationPval] = para.test(sc(:,1:para.nPC), ...
                    covariate,'type','spearman','rows','pairwise');
            else
                report.pcCovAssociationPval = NaN(para.nPC, nvar);
                for i = 1:para.nPC
                    for j = 1:nvar               
                        report.pcCovAssociationPval(i,j) = para.test(sc(:,i),covariate(:,j),'off');
                    end
                end
            end
        end
        
        function [fig, cbarh] = plot(report, covIdx, varargin)
            %report: from test function
            %covIdx: index to report.covariate or a string (name of the
            %   covariate)
            para.fig = [];
            para.nPC = 10;
            para.legendAll = false;
            para.colormap = genColorMap([1, 0.8, 0.3; 0.3, 0.3, 0.8],10);
            para = assignpara(para, varargin{:});            
            para.nPC = min(para.nPC, report.para.nPC);
            
            cbarh = NaN;
            if ischar(covIdx)
                [~, i] = ismember(covIdx, report.para.covLabel);
                if i == 0
                    error('cannot find %s', covIdx);
                end
                covIdx = i;
            end            
            if ~isempty(para.fig)
                set(0, 'CurrentFigure', para.fig);
                fig = para.fig;
            else
                fig = gcf;
            end
            
            if strcmp(func2str(report.para.test), 'corr')
                [~, sortDataIdx] = sort(report.covariate(:,covIdx),'descend');
            end
            
            [plotRow, plotCol] = numSubplot(para.nPC-1);
            [~, sortDataIdx] = sort(report.covariate(:, covIdx));
            for i = 1:para.nPC-1
                subplot(plotRow, plotCol, i);
                if strcmp(func2str(report.para.test), 'corr')                    
                    scatter(report.pcProjection(sortDataIdx,i), report.pcProjection(sortDataIdx,i+1), ...
                        50, report.covariate(sortDataIdx,covIdx), 'fill');
                    colormap(para.colormap);
                    if i == 1
                        apos = get(gca, 'position');
                    end
                    if i == para.nPC-1
                        cbarh = colorbar;
                        cpos = get(cbarh, 'position');
                        set(cbarh, 'position', [0.93, 0.4, cpos(3), cpos(4)]);
                        curapos = get(gca, 'position');
                        set(gca, 'position', [curapos(1), curapos(2), apos(3), curapos(3)]);
                    end
                else
                    h = gscatter(report.pcProjection(sortDataIdx,i), ...
                        report.pcProjection(sortDataIdx,i+1), ...
                        report.covariate(sortDataIdx,covIdx));
                    distinctColorMarker(h);
                    if ~para.legendAll && i ~= 1%para.nPC-1
                        legend off; 
                    end
                end
                xlim([min(report.pcProjection(:,i))-0.1, max(report.pcProjection(:,i)+0.1)]);
                ylim([min(report.pcProjection(:,i+1))-0.1, max(report.pcProjection(:,i+1)+0.1)]);
                xlabel(sprintf('PC%d, %0.1f%%',i,report.pcVarExplained(i)*100),'fontsize',10);
                ylabel(sprintf('PC%d, %0.1f%%',i+1,report.pcVarExplained(i+1)*100),'fontsize',10);
                title(sprintf('%s pval w/ x,y:\n%0.1e, %0.1e',func2str(report.para.test), ...
                    report.pcCovAssociationPval(i,covIdx), ...
                    report.pcCovAssociationPval(i+1,covIdx)), 'fontsize',10);                                
            end
            suptitle(sprintf('Sample labeled by %s', strrep(report.para.covLabel{covIdx}, '_', ' ')));
            lh = findobj(fig, 'type', 'axes', 'tag', 'legend');
            uistack(lh, 'top');
        end
        
        function fig = linearplot(report, covIdx, varargin)
            %report: from test function
            %covIdx: index to report.covariate or a string (name of the
            %   covariate)
            para.fig = [];
            para.nPC = 9;                       
            para = assignpara(para, varargin{:});            
            para.nPC = min(para.nPC, report.para.nPC);            
            
            if ischar(covIdx)
                [~, i] = ismember(covIdx, report.para.covLabel);
                if i == 0
                    error('cannot find %s', covIdx);
                end
                covIdx = i;
            end            
            if ~isempty(para.fig)
                set(0, 'CurrentFigure', para.fig);
                fig = para.fig;
            else
                fig = gcf;
            end
                        
            [plotRow, plotCol] = numSubplot(para.nPC);
            [~, sortDataIdx] = sort(report.covariate(:, covIdx));
            for i = 1:para.nPC
                subplot(plotRow, plotCol, i);
                plot(report.pcProjection(sortDataIdx,i), report.covariate(sortDataIdx,covIdx),'.');                
                xlim([min(report.pcProjection(:,i))-0.1, max(report.pcProjection(:,i))+0.1]);
                ylim([min(report.covariate(sortDataIdx,covIdx))-0.1, max(report.covariate(sortDataIdx,covIdx))+0.1]);
                xlabel(sprintf('PC%d, %0.1f%%',i,report.pcVarExplained(i)*100),'fontsize',10);
                ylabel(sprintf('covariate'),'fontsize',10);
                title(sprintf('%s pval: %0.1e',func2str(report.para.test), ...
                    report.pcCovAssociationPval(i,covIdx)), 'fontsize',10);                                
            end
            suptitle(sprintf('Sample labeled by %s', strrep(report.para.covLabel{covIdx}, '_', ' ')));
        end
        
        function fig = boxplot(report, covIdx, varargin)
            %report: from test function
            %covIdx: index to report.covariate or a string (name of the
            %   covariate)
            para.fig = [];
            para.nPC = 9;
            para.legendAll = false;
            para.showBatchLabel = true;
            para = assignpara(para, varargin{:});            
            para.nPC = min(para.nPC, report.para.nPC);
                        
            if ischar(covIdx)
                [~, i] = ismember(covIdx, report.para.covLabel);
                if i == 0
                    error('cannot find %s', covIdx);
                end
                covIdx = i;
            end            
            if ~isempty(para.fig)
                set(0, 'CurrentFigure', para.fig);
                fig = para.fig;
            else
                fig = gcf;
            end
            
            [plotRow, plotCol] = numSubplot(para.nPC);
            [~, covSortIdx] = sort(report.covariate(:, covIdx));
            for i = 1:para.nPC
                subplot(plotRow, plotCol, i);
                boxplotWithData(report.pcProjection(covSortIdx,i), report.covariate(covSortIdx,covIdx));
%                 boxplot(report.pcProjection(covSortIdx,i), report.covariate(covSortIdx,covIdx));                                                
                ylabel(sprintf('PC%d, %0.1f%%',i,report.pcVarExplained(i)*100),'fontsize',12);
                title(sprintf('%s pval: %0.1e',func2str(report.para.test), ...
                    report.pcCovAssociationPval(i,covIdx)), 'fontsize',12);                                
                if ~para.showBatchLabel
                    set(gca, 'xticklabel', {' '});
                end
            end
            suptitle(sprintf('Distribution of PC by %s', report.para.covLabel{covIdx}));
        end
        
        function table = makeTable(report)
                        
            [nPC, nCov] = size(report.pcCovAssociationPval);
            n = nPC * nCov;

            table = cell(n+1, 4);
            table(1,:) = {'PC','variance explained (%)', 'covariate', ''};
            table{1,4} = sprintf('%s pval',func2str(report.para.test));            
            table(2:end,1) = numarray2strarray( reshape(repmat((1:nPC)',1,nCov), n, 1));
            table(2:end,2) = numarray2strarray( reshape(repmat(report.pcVarExplained*100,1,nCov),n,1));
            table(2:end,3) = reshape(repmat(report.para.covLabel(:)',nPC,1),n,1);
            table(2:end,4) = numarray2strarray( reshape( report.pcCovAssociationPval, n, 1));
            if length(unique(table(2:end,3))) == 1
                table(:,3) = [];
            end
        end
    end
end
            

    