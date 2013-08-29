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
function [pvalue,score,cutoff,num_annotated,num_non_annotated] = ...
    GSEA_old(r, annotations, score_func, P, num_permutations)

    %r: #gene x 1
    %annotations: #gene x #annotation
    %

    if(nargin < 4)
        P = 1;
    end
    if(nargin < 5)
        num_permutations = 1000;
    end

    [numgene numannotation] = size(annotations);
    
    [sortedR, rOrder] = sort(r,1,'descend');
    sortedAnnotations = annotations(rOrder, :); %#gene x #annotation
    member_scores = feval(score_func, sortedR); %#gene x 1

%     NR = sum(abs(member_scores(sortedAnnotations==1)).^P);    
    %use for loop instead of bsxfun or repmat, because repmat require
    %memory and bsxfun might be problematic with sparse matrix (at least on
    %2010b version)
    NR = full(sum(bsxfun(@times, abs(member_scores), sortedAnnotations==1) .^P, 1)); %1 x #annotation

    [score,num_annotated,num_non_annotated] = computeScore(member_scores, sortedAnnotations, P, NR);
    cutoff(1,:) = sortedR(num_annotated(1,:) + num_non_annotated(1,:));
    cutoff(2,:) = sortedR(length(sortedR) - num_annotated(2,:) - num_non_annotated(2,:)+1);

    distribution1 = zeros(num_permutations, numannotation);
    distribution2 = zeros(num_permutations, numannotation);
    for perm = 1 : num_permutations
        random_annotation = annotations(randperm(numgene), :);
%         NRrand = sum(abs(member_scores(random_annotation==1)).^P);
        NRrand = full(sum(bsxfun(@times, abs(member_scores), random_annotation==1).^P, 1));
        sc = computeScore(member_scores, random_annotation, P, NRrand);
        distribution1(perm,:) = sc(1,:);
        distribution2(perm,:) = sc(2,:);
    end

    pvalue = [ sum(distribution1 > repmat(score(1,:), num_permutations, 1) ); ...
               sum(distribution2 > repmat(score(2,:), num_permutations, 1) )] ./ num_permutations;
end


function [score,num_annotated,num_non_annotated] = computeScore(member_scores, annotated, P, NR)
    %member_scores: #gene x 1
    %annotated: #gene x #annotation

    miss_penalty = 1 ./ sum(annotated == 0,1); %1 x #annotation


    score = zeros(2, size(annotated,2));
    num_annotated = NaN(2, size(annotated,2));
    num_non_annotated = NaN(2, size(annotated,2));
    
    tmpscore = bsxfun(@rdivide, member_scores.^P, NR); %#gene x annotation
    tmpscore(annotated ~= 1) = 0; %just in case there are values other than 0 and 1
    tmp2 = repmat(-miss_penalty, size(tmpscore,1), 1);
    tmpscore(annotated == 0) = tmp2(annotated==0);
    clear tmp2
    
    [score(1,:), i] = max(cumsum(tmpscore,1), [], 1);
    for j = 1:size(annotated,2)
        num_annotated(1,j) = sum(annotated(1:i(j),j)==1, 1); % 1 x #annotation
        num_non_annotated(1,j) = sum(annotated(1:i(j),j)==0, 1);
    end
    
    score(1,:) = max(score(1,:),0);
    num_annotated(1, score(1,:)==0) = 0;
    num_non_annotated(1, score(1,:)==0) = sum(annotated(:,score(1,:)==0),1);
    
    tmpscore = tmpscore(size(tmpscore,1):-1:1,:);
    rannotated = annotated(size(tmpscore,1):-1:1, :);
    [score(2,:), i] = max(cumsum(tmpscore,1), [], 1);
    for j = 1:size(annotated,2)
        num_annotated(2,j) = sum(rannotated(1:i(j),j)==1, 1); % 1 x #annotation
        num_non_annotated(2,j) = sum(rannotated(1:i(j),j)==0, 1);
    end
    score(2,:) = max(score(2,:), 0);
    num_annotated(2, score(2,:)==0) = 0;
    num_non_annotated(2, score(2,:)==0) = sum(rannotated(:,score(2,:)==0),1);
    clear rannotated tmpscore
        
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