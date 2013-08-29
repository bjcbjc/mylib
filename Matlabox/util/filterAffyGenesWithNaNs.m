function [filtered_data, kept_genes, kept_probes,removed_genes] = filterAffyGenesWithNaNs(genes,data, correlation_threshold)
    % This script will filter expression results, starting with the raw probe
    % level, and ending up with gene level
    keep_probes = ones(size(genes));
    kept_probes = false(size(data,1),1);
    if (nargin < 3)
        correlation_threshold = 0.75;
    end
    numsamples=size(data,2);
    
    rand('state', 1);

    uniquegenes=unique(genes);
    filtered_data=zeros(length(uniquegenes),numsamples);
    keep_genes=zeros(length(uniquegenes),1);
    for i=1:length(uniquegenes)
        gene=uniquegenes(i);
        transcripts = (genes == gene);
        keep_gene = 0;
        gene_data = [];
        gene_expression = data(transcripts, :);
        if(sum(transcripts) == 1)
            keep_gene = 1;
            gene_data = data(transcripts, :);
            kept_probes(transcripts) = true;
        elseif (sum(any(isnan(gene_expression),2)) == sum(transcripts)) % all transcripts have NaNs
            [~, b] = min(sum(isnan(gene_expression), 2));
            [~, d] = max(nanmean(gene_expression, 2));
            minCorr = min(corr(gene_expression',nanmean(gene_expression)', 'rows', 'complete'));
            if (b == d && min(sum(isnan(gene_expression), 2)) ~= max(sum(isnan(gene_expression), 2)))
                gene_data = gene_expression(b, :);
                keep_gene = 1;
                tr = find(transcripts==1);
                kept_probes(tr(b)) = true;
            elseif (minCorr > correlation_threshold)
                keep_gene = 1;
                gene_data = mean(gene_expression);
                kept_probes(transcripts) = true;
            else % pick the one with the higher values
                keep_gene = 1;
                [~,mi] = max(nanmean(gene_expression));
                gene_data = gene_expression(mi,:);
                ps = find(transcripts);
                transcripts(ps(3-mi)) = false;
                kept_probes(transcripts) = true;
            end
        elseif (sum(any(isnan(gene_expression),2)) == sum(transcripts) - 1) % one transcript without NaNs - use it
            gene_data = gene_expression(~any(isnan(gene_expression),2), :);
            tr = find(transcripts==1);
            kept_probes(tr(~any(isnan(gene_expression),2))) = true;
            keep_gene = 1;
        else
            gene_expression = gene_expression(~any(isnan(gene_expression),2), :);
            tr = find(transcripts==1);
            tr = tr(~any(isnan(gene_expression),2));
            %minCorr = min(min(corr(gene_expression')));
            minCorr = min(corr(gene_expression',nanmean(gene_expression)'));

            if (minCorr > correlation_threshold)
                keep_gene = 1;
                gene_data = mean(gene_expression);
                kept_probes(tr) = true;
            else
                [without_low,keptprobes] = without_low_probes(gene_expression, correlation_threshold);
                kept_probes(tr(keptprobes)) = true;
                %if (size(without_low,1) == 1)
                keep_gene = 1;
                gene_data = without_low;                
                %elseif (size(without_low,1) > 1)
                %    [gene_data, keep_gene] = filterNoise(without_low, correlation_threshold);
                %end
%                 if (keep_gene == 0)
%                     spearman = corr(gene_expression', 'type', 'Spearman');
%                     spearman_correlated = spearman > correlation_threshold;
%                     if (max(sum(spearman_correlated)) > 1)
%                         maxcorrprobe=find(sum(spearman_correlated)==max(sum(spearman_correlated)));
%                         gene_expression_correlated = gene_expression(spearman_correlated(maxcorrprobe(1),:),:);
%                         % After looking at the genes where only one is left, I
%                         % don't like them
%                         if (size(gene_expression_correlated, 1) > 1)
%                             [gene_data, keep_gene] = filterNoise(gene_expression_correlated, correlation_threshold);
%                         else
%                             %selectedprobes(transcripts(~spearman_correlated(maxcorrprobe(1),:)))=[];
%                             keep_probes(transcripts(~spearman_correlated(maxcorrprobe(1),:))) = 0;
%                         end
%                     else
%                         [gene_data, keep_gene] = filterNoise(gene_expression, correlation_threshold);
%                         if (keep_gene == 0)
%                             spearman = corr(gene_expression', 'type', 'Spearman');
%                             spearman_correlated = spearman > correlation_threshold;
%                             if (max(sum(spearman_correlated)) > 1)
%                                 maxcorrprobe=find(sum(spearman_correlated)==max(sum(spearman_correlated)));
%                                 gene_expression_correlated = gene_expression(spearman_correlated(maxcorrprobe(1),:),:);
%                                 % After looking at the genes where only one is left, I
%                                 % don't like them
%                                 if (size(gene_expression_correlated, 1) > 1)
%                                     [gene_data, keep_gene] = filterNoise(gene_expression_correlated, correlation_threshold);
%                                 end
%                             end
%                         end
%                     end
%                 end
            end
        end
        if (keep_gene)
            % Skip Mean normalization
            %valid = not(isnan(gene_data));
            %gene_data = gene_data - mean(gene_data(:, valid));
            filtered_data(i,:) = gene_data;
            %[filtered_data; gene_data];
            keep_genes(i)=1;
        else
            keep_genes(i)=0;
            %selectedprobes(transcripts)=[];
            keep_probes(transcripts) = 0;
        end
    end
    kept_genes=uniquegenes(keep_genes==1);
    removed_genes=uniquegenes(keep_genes==0);
    filtered_data=filtered_data(keep_genes==1,:);
end

function [filteredExpression, keepGene] = filterNoise(geneExpression, correlation_threshold)        
    filteredExpression = geneExpression;
    keepGene = 0;
    [geneExpressionFiltered, samplesToKeep] = filter_noisy_samples(geneExpression, correlation_threshold);
    %minCorr = min(min(corr(geneExpressionFiltered(:, samplesToKeep)')));
    minCorr=min(corr(geneExpressionFiltered(:, samplesToKeep)',nanmean(geneExpressionFiltered(:, samplesToKeep))'));
    if (minCorr > correlation_threshold)
        keepGene = 1;
        filteredExpression = mean(geneExpressionFiltered);        
    end
    normalized = (geneExpression - repmat(mean(geneExpression, 2), 1, size(geneExpression, 2))) ./ repmat(std(geneExpression, 0, 2), 1, size(geneExpression, 2));
    % I'm not using abs > 3, because of this reason:
    % When I have a probe that's very low, it makes biological sense that it is
    % the outlier.
    % When I have a probe that's very high, I don't know if the biological
    % "truth" is that there was a high expression value in these samples, and
    % the other probe didn't respond, or that there was a low expression value
    % in these samples, and that the high values are the outliers.
    probes_with_outliers = any(normalized(:, not(samplesToKeep))  < -3, 2);
    if (sum(probes_with_outliers) < size(geneExpression, 1)/2)
        geneExpression = geneExpression(not(probes_with_outliers), :);
        %minCorr = min(min(corr(geneExpression')));
        minCorr=min(corr(geneExpression',nanmean(geneExpression)'));
        if (minCorr > correlation_threshold)
            keepGene = 1;
            filteredExpression = mean(geneExpression);            
        end
    end
end

function [gene_expression, to_keep] = filter_noisy_samples(gene_expression, correlation_threshold, max_fraction_to_remove)
    if (nargin < 2)
        correlation_threshold = 0.75;
    end
    if (nargin < 3)
        max_fraction_to_remove = 1/15;
    end

    max_num_to_remove = floor(length(gene_expression(1, :)) * max_fraction_to_remove);

    to_keep = true(size(gene_expression(1, :)));
    %if (min(min(corr(gene_expression'))) < correlation_threshold)
    if (min(corr(gene_expression',nanmean(gene_expression)')) < correlation_threshold)
        % deviations = sort(std(gene_expression), 'descend');
        num_cols = length(gene_expression(1, :));
        meanStdNormalized = (gene_expression - repmat(mean(gene_expression, 2), 1, num_cols)) ./ repmat(std(gene_expression, 0, 2), 1, num_cols);
        deviations = sort(std(meanStdNormalized), 'descend');
        for deviation = deviations(1:(max_num_to_remove + 1))
        %for deviation = deviations(1:end-2)
            to_remove = std(meanStdNormalized) >= deviation;
            to_keep = not(to_remove);
            %if (min(min(corr(gene_expression(:, to_keep)'))) > correlation_threshold)
            if (min(corr(gene_expression(:, to_keep)',nanmean(gene_expression(:, to_keep))')) > correlation_threshold)
                break
            end
        end
        if (sum(to_remove) <= max_num_to_remove)
            gene_expression(:, to_remove) = NaN;
        else
            to_keep = true(size(gene_expression(1, :)));
        end
    end
end

function [gene_expression, kept_probes] = without_low_probes(gene_expression, min_corr_th)
% Filters out probes whose mean is at least N stdevs lower (default N =2)
% than the mean of the entire gene expression matrix.
% gene_expression : Matrix of gene expression where columns are samples,
%                   and each row is one probe. All the probes are supposed
%                   to belong to the same gene.
c = min(corr(gene_expression',nanmean(gene_expression)'));
kept_probes = 1:size(gene_expression,1);
while c<min_corr_th && size(gene_expression,1)>1
    if (size(gene_expression,1)>2)
        c = corr(gene_expression');
        m=min(c);
        [~,mi]=min(m);
        kept_probes = kept_probes((1:size(gene_expression,1))~=mi);
        gene_expression = gene_expression( (1:size(gene_expression,1))~=mi ,:);
        c = min(corr(gene_expression',nanmean(gene_expression)'));
    else        
        [~,mi]=min(mean(gene_expression, 2));
        kept_probes = kept_probes((1:size(gene_expression,1))~=mi);
        gene_expression = gene_expression( (1:size(gene_expression,1))~=mi ,:);
        if (size(gene_expression,1)>1)
            c = min(corr(gene_expression',nanmean(gene_expression)'));
        end
    end
end
if (size(gene_expression,1)>1)
    gene_expression = mean(gene_expression);
end
end