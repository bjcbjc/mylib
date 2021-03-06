%rankedGSEA Run GSEA algorithm using rankings
%
% function [pvalue,score,cutoff,num_annotated,num_non_annotated] = GSEA(r,annotations, score_func, P, num_permutations)
%
% The two first input parameters, r and annotations, must be row
% vectors of the same length. r provides the score of each gene;
% annotations is 1 for each annotated gene, 0 for each non-annotated one.
%
% score_func - the function used to translate the r values to a score for each 
% annotated member in the set. You can write different functions and pass them
% as an argument to test different scoring methods. The function receives a row 
% vector as an argument, which it can assume are the r values sorted in descending
% order; and returns a row vector of the same length, comparing the score values to
% use in the GSEA algorithm. To use the absolute value of the r value as the score 
% (the method described in the GSEA paper, and used in our currently-implemented 
% java version) use 'abs' as the value of score_func.
%
% P (optional argument) - raise the results of score_func to the power of P. 
% P usually should be between 0 and 1. Default is 1. Note that if P is 0, it 
% makes no difference what score_func is.
%
% num_permutations (optional) - default is 1000
%
% Each of the return arguments is a vector of length 2; the first value is
% for the result of the forward search, the second for the result of the
% reverse search.
%
% modified by BJ: take multiple annotations
%
%
function [pvalue, qvalue, score,cutoff,num_annotated,num_non_annotated] = ...
    GSEA(r, annotations, score_func, P, num_permutations, varargin)

    %r: #gene x 1
    %annotations: #gene x #annotation
    %
    para.quiet = false;
    para.direction = [1, 1]; % large-value, small-value enrichment
    
    if(nargin < 4)
        P = 1;
    end
    if(nargin < 5)
        num_permutations = 1000;
    end

    para = assignpara(para, varargin{:});
    
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
    
    
    %this is a costly line; calculate outsie
    tmpscore = member_scores .^ P;
    %tmpscore = bsxfun(@rdivide, member_scores.^P, NR); %#gene x annotation
    
    member_scores = abs(member_scores);
    
%     NR = sum(abs(member_scores(sortedAnnotations==1)).^P);    
    %use for loop instead of bsxfun or repmat, because repmat require
    %memory and bsxfun might be problematic with sparse matrix (at least on
    %2010b version)
    %NR = full(sum(bsxfun(@times, abs(member_scores), sortedAnnotations==1) .^P, 1)); %1 x #annotation
    NR = sum(bsxfun(@times, member_scores, sortedAnnotations) .^P, 1); %1 x #annotation

    %this will not change during permutation
    miss_penalty = 1 ./ sum(~annotations, 1); %1 x #annotation
    miss_penalty = repmat(-miss_penalty, numgene, 1);    %gene x annotation
    
    [score,num_annotated,num_non_annotated] = computeScore(tmpscore, sortedAnnotations, miss_penalty, NR, para);
    cutoff = NaN(2, numannotation);
    if para.direction(1)
        cutoff(1,:) = sortedR(num_annotated(1,:) + num_non_annotated(1,:));
    end
    if para.direction(2)
        cutoff(2,:) = sortedR(numgene - num_annotated(2,:) - num_non_annotated(2,:)+1);
    end

    distribution1 = zeros(num_permutations, numannotation);
    distribution2 = zeros(num_permutations, numannotation);
    
    if num_permutations < numannotation 
        for perm = 1 : num_permutations
            if mod(perm, round(num_permutations/10)) == 0 && ~para.quiet
                fprintf('%d permutations done\n', perm);
            end
            random_annotation = annotations(randperm(numgene), :);
    %         NRrand = sum(abs(member_scores(random_annotation==1)).^P);
            NRrand = sum(bsxfun(@times, member_scores, random_annotation).^P, 1);
            sc = computeScore(tmpscore, random_annotation, miss_penalty, NRrand, para);
            distribution1(perm,:) = sc(1,:);
            distribution2(perm,:) = sc(2,:);
        end
    else
        permIdx = NaN(numgene, num_permutations);
        for perm = 1:num_permutations
            permIdx(:, perm) = randperm(numgene);
        end        
        for anti = 1:numannotation            
            random_annotation = reshape(annotations(permIdx, anti), numgene, num_permutations);    
            miss_penalty_perm = repmat(miss_penalty(:, anti), 1, num_permutations);
            NRrand = sum(bsxfun(@times, member_scores, random_annotation).^P, 1);
            sc = computeScore(tmpscore, random_annotation, miss_penalty_perm, NRrand, para);
            distribution1(:,anti) = sc(1,:);
            distribution2(:,anti) = sc(2,:);
        end
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
    rmi = find(para.direction == 0);
    pvalue(rmi,:) = [];
    qvalue(rmi,:) = [];
    score(rmi, :) = [];
    cutoff(rmi,:) = [];
    num_annotated(rmi,:) = [];
    num_non_annotated(rmi,:) = [];
%     pvalue = [ sum(distribution1 > repmat(score(1,:), num_permutations, 1) ); ...
%                sum(distribution2 > repmat(score(2,:), num_permutations, 1) )] ...
%                ./ num_permutations;
end


function [score,num_annotated,num_non_annotated] = computeScore(tmpscore, annotated, defaultscore, NR, para)
    %member_scores: #gene x 1
    %annotated: #gene x #annotation

    %miss_penalty = 1 ./ sum(annotated == 0,1); %1 x #annotation

    [ngene, nant] = size(annotated);
    
    score = zeros(2, nant);
    num_annotated = NaN(2, nant);
    num_non_annotated = NaN(2, nant);
    
    tmp2 = bsxfun(@rdivide, tmpscore, NR);
    %tmpscore = bsxfun(@rdivide, member_scores.^P, NR); %#gene x annotation
    %tmpscore(~annotated) = 0; %just in case there are values other than 0 and 1
    
    defaultscore(annotated) = tmp2(annotated);
    clear tmp2
    
    if para.direction(1)
        [score(1,:), i] = max(cumsum(defaultscore,1), [], 1);
        score(1,:) = max(score(1,:),0);

        if nargout > 1
            I = repmat( (1 : ngene)', 1, nant);
            idx = bsxfun(@le, I, i);
            num_annotated(1,:) = sum( annotated .* idx, 1);
            %     cp = annotated(idx);
            %     annotated(idx) = 0;
            %     num_annotated(1,:) = sum(annotated, 1);
            num_non_annotated(1,:) = i - num_annotated(1,:);

            %     annotated(idx) = cp; %recover annotated

            %     for j = 1:nant
            %         num_annotated(1,j) = sum(annotated(1:i(j),j)==1, 1); % 1 x #annotation
            %         num_non_annotated(1,j) = sum(annotated(1:i(j),j)==0, 1);
            %     end

            num_annotated(1, score(1,:)==0) = 0;
            num_non_annotated(1, score(1,:)==0) = sum(annotated(:,score(1,:)==0),1);
        end
    end
       
    if para.direction(2)
        defaultscore = defaultscore(ngene:-1:1, :);
        rannotated = annotated(ngene:-1:1, :);
        [score(2,:), i] = max(cumsum(defaultscore,1), [], 1);
        score(2,:) = max(score(2,:), 0);

        if nargout > 1
            idx = bsxfun(@le, I, i);
            num_annotated(2,:) = sum( rannotated .* idx, 1);
            %     cp = rannotated(idx);
            %     rannotated(idx) = 0;
            %     num_annotated(2,:) = sum(rannotated, 1);
            num_non_annotated(2,:) = i - num_annotated(2,:);

            %     rannotated(idx) = cp;

            %     for j = 1:nant
            %         num_annotated(2,j) = sum(rannotated(1:i(j),j)==1, 1); % 1 x #annotation
            %         num_non_annotated(2,j) = sum(rannotated(1:i(j),j)==0, 1);
            %     end        
            num_annotated(2, score(2,:)==0) = 0;
            num_non_annotated(2, score(2,:)==0) = sum(rannotated(:,score(2,:)==0),1);
        end
    end
end

function input = returnSelf(input)
end

function rankings = NMranking(r)
% function NMranking
%
% This is a function intended to be used as the score_func argument to
% GSEA. It scores each member by its rank within the positive or negative
% r values. If there are N members with positive r values and M members with
% negative values, then the scores of the positive members will be 1 thru
% N, with the highest member having score N; and the scores of the
% non-positive members will be 1 thru M, with the lowest member having score M
    N = length(find(r > 0));
    M = length(r) - N;
    rankings = [N : -1 : 1 1 : M];
end