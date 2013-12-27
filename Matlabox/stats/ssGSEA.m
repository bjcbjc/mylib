function [score, pvalue, qvalue] = ssGSEA(r, annotations, score_func, P, num_permutations, reverseSearch)
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
% Each of the return arguments is a vector of length 2 (if reverse search is on); 
%   the first value is for the result of the forward search, the second for the result of the
%   reverse search. 
%

    %r: #gene x 1
    %annotations: #gene x #annotation
    %

    if nargin < 3, score_func = @abs; end
    if nargin < 4, P = 1; end
    if nargin < 5, num_permutations = 0; end
    if nargin < 6, reverseSearch = true; end

    [numgene, numannotation] = size(annotations);
    
    if issparse(annotations)
        annotations = full(annotations);
    end
    if ~islogical(annotations)
        annotations = annotations ~= 0;
    end
    
    [sortedR, rOrder] = sort(r,1,'descend');
%     [sortedR, rOrder] = sort(r,1,'ascend');
    
    sortedAnnotations = annotations(rOrder, :); %#gene x #annotation
    member_scores = feval(score_func, sortedR); %#gene x 1
    
    
    %this is a costly line; calculate outside
    tmpscore = member_scores .^ P;
    %tmpscore = bsxfun(@rdivide, member_scores.^P, NR); %#gene x annotation
        
    NR = sum(bsxfun(@times, member_scores, sortedAnnotations) .^P, 1); %1 x #annotation

    %this will not change during permutation
    miss_penalty = 1 ./ sum(~annotations, 1); %1 x #annotation
    miss_penalty = repmat(-miss_penalty, numgene, 1);    
    
    score = computeScore(tmpscore, sortedAnnotations, miss_penalty, NR, reverseSearch);
    
    pvalue = [];
    qvalue = [];
    if num_permutations > 0
        distribution1 = zeros(num_permutations, numannotation, 1+reverseSearch);
        
        for perm = 1 : num_permutations
%             if mod(perm, round(num_permutations/10)) == 0
%                 fprintf('%d permutations done\n', perm);
%             end
            random_annotation = annotations(randperm(numgene), :);
            NRrand = sum(bsxfun(@times, member_scores, random_annotation).^P, 1);
            sc = computeScore(tmpscore, random_annotation, miss_penalty, NRrand, reverseSearch);
            for searchIdx = 1:reverseSearch+1
                distribution1(perm,:,searchIdx) = sc(searchIdx,:);
            end
        end
        %normalization of score; BJ, 11/09/2012
        meanscore = mean(distribution1, 1); % 1 x #annotation x #search (forward+reverse?)
        if ndims(meanscore) == 3
            meanscore = squeeze(meanscore)'; % #search x #annotation
        end
        normscore = score ./ meanscore; % #search x #annotation
        rcount = zeros(1+reverseSearch, numannotation);
        psratio = zeros(1+reverseSearch, numannotation);
        for searchIdx = 1:reverseSearch+1
            distribution1(:,:,searchIdx) = bsxfun(@rdivide, distribution1(:,:,searchIdx), meanscore(searchIdx,:));
            nulldist = distribution1(:,:,searchIdx);            
            rcount(searchIdx, :) = sum( bsxfun(@ge, nulldist(:), normscore(searchIdx,:)), 1); 
            psratio(searchIdx, :) = mean(bsxfun(@ge, normscore(searchIdx,:)', normscore(searchIdx,:)), 1);
        end
        pvalue = rcount ./ num_permutations ./ numannotation;        
        qvalue = pvalue ./ psratio;
        %     pvalue = [ sum(distribution1 > repmat(score(1,:), num_permutations, 1) ); ...
        %                sum(distribution2 > repmat(score(2,:), num_permutations, 1) )] ...
        %                ./ num_permutations;    
    end
end


function score = computeScore(tmpscore, annotated, defaultscore, NR, revsearch)
    %tmpscores: #gene x #annotation
    %annotated: #gene x #annotation
    %defaultscore: #gene x #annotation, -P_miss
    %NR: 1 x #annotation, normalizing factor for hits

    [ngene, nant] = size(annotated);
    
    score = zeros(1+revsearch, nant);
    
    tmp2 = bsxfun(@rdivide, tmpscore, NR);
    
    defaultscore(annotated) = tmp2(annotated);
    clear tmp2
    
    % different from GSEA which takes max of deviation as enrichment score,
    % ssGSEA takes the sum 
    score(1,:) = sum(cumsum(defaultscore, 1), 1);
    %original GSEA
%     score(1,:) = max(cumsum(defaultscore, 1), [], 1);
    if revsearch
        defaultscore = defaultscore(ngene:-1:1, :);
        score(2,:) = sum(cumsum(defaultscore, 1), 1);
    end
end


