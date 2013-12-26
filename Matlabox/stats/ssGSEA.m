function [pvalue, qvalue, score] = ssGSEA(r, annotations, score_func, P, num_permutations)
% ssGSEA: single sample GSEA
% use ranking of data in a single sample to calculate enrichment in gene
% sets
%
% modified from GSEA, originally coded by Dylan, modified by BJ
%
% input:
%   r: #genes x 1; r provides the score of each gene 
%   annotations: #genes x #annotation; 1 for each annotated gene, 0 for each non-annotated one.
%   score_func: the function used to translate the r values to a score for each 
%       annotated member in the set. default = @abs
%   P: raise the results of score_func to the power of P; default = 0.25
%       P usually should be between 0 and 1. Default is 1. Note that if P is 0, it 
%       makes no difference what score_func is.
%   num_permutations: default = 1000
%
% Each of the return arguments is a vector of length 2; the first value is
% for the result of the forward search, the second for the result of the
% reverse search.
%

    %r: #gene x 1
    %annotations: #gene x #annotation
    %

    if nargin < 3, score_func = @abs; end
    if nargin < 4, P = 1; end
    if nargin < 5, num_permutations = 1000; end

    [numgene, numannotation] = size(annotations);
    
    if issparse(annotations)
        annotations = full(annotations);
    end
    if ~islogical(annotations)
        annotations = annotations ~= 0;
    end
    
    [sortedR, rOrder] = sort(r,1,'descend');
    sortedAnnotations = annotations(rOrder, :); %#gene x #annotation
    member_scores = feval(score_func, sortedR); %#gene x 1
    
    
    %this is a costly line; calculate outside
    tmpscore = member_scores .^ P;
    %tmpscore = bsxfun(@rdivide, member_scores.^P, NR); %#gene x annotation
        
    NR = sum(bsxfun(@times, member_scores, sortedAnnotations) .^P, 1); %1 x #annotation

    %this will not change during permutation
    miss_penalty = 1 ./ sum(~annotations, 1); %1 x #annotation
    miss_penalty = repmat(-miss_penalty, numgene, 1);    
    
    score = computeScore(tmpscore, sortedAnnotations, miss_penalty, NR);
    
    distribution1 = zeros(num_permutations, numannotation);
    distribution2 = zeros(num_permutations, numannotation);
    for perm = 1 : num_permutations
        if mod(perm, round(num_permutations/10)) == 0
            fprintf('%d permutations done\n', perm);
        end
        random_annotation = annotations(randperm(numgene), :);
        NRrand = sum(bsxfun(@times, member_scores, random_annotation).^P, 1);
        sc = computeScore(tmpscore, random_annotation, miss_penalty, NRrand);
        distribution1(perm,:) = sc(1,:);
        distribution2(perm,:) = sc(2,:);
    end
    %normalization of score; BJ, 11/09/2012
    meanscore = [mean(distribution1, 1); mean(distribution2, 1)];
    distribution1 = bsxfun(@rdivide, distribution1, meanscore(1,:));
    distribution2 = bsxfun(@rdivide, distribution2, meanscore(2,:));
    normscore = score ./ meanscore;
        
    rcount = [sum( bsxfun(@ge, distribution1(:), normscore(1,:)), 1); ...
        sum(bsxfun(@ge, distribution2(:), normscore(2,:)), 1)] ;
    pvalue = rcount ./ num_permutations ./ numannotation;
    psratio = [mean(bsxfun(@ge, normscore(1,:)', normscore(1,:)), 1); ...
        mean(bsxfun(@ge, normscore(2,:)', normscore(2,:)), 1)] ;
    qvalue = pvalue ./ psratio;
%     pvalue = [ sum(distribution1 > repmat(score(1,:), num_permutations, 1) ); ...
%                sum(distribution2 > repmat(score(2,:), num_permutations, 1) )] ...
%                ./ num_permutations;
end


function score = computeScore(tmpscore, annotated, defaultscore, NR)
    %tmpscores: #gene x #annotation
    %annotated: #gene x #annotation
    %defaultscore: #gene x #annotation, -P_miss
    %NR: 1 x #annotation, normalizing factor for hits

    [ngene, nant] = size(annotated);
    
    score = zeros(2, nant);
    
    tmp2 = bsxfun(@rdivide, tmpscore, NR);
    
    defaultscore(annotated) = tmp2(annotated);
    clear tmp2
    
    % different from GSEA which takes max of deviation as enrichment score,
    % ssGSEA takes the sum 
    score(1,:) = sum(cumsum(defaultscore, 1), 1);
    %original GSEA
    %[score(1,:), i] = max(cumsum(defaultscore,1), [], 1);
    %score(1,:) = max(score(1,:),0);
        
    defaultscore = defaultscore(ngene:-1:1, :);
    score(2,:) = sum(cumsum(defaultscore, 1), 1);        
end


