function [precision, recall, fmeasure] = retrieval(answer, guess, correct_zerofeature, fbeta)
%answer and guess are logicals
%if answer and guess are matrices, the dimension should be #feature x
%#response
%
%nonan: if sum(answer)==sum(guess)==0, 0/0 = NaN (precision or recall). 
%   correct_zerofeature = true: correct this to 1 (because correctly guess no feaure)
%   correct_zerofeature = false: leave it as NaN
%
%fbeta: used for fmeasure = (1+fbeta^2)*precision*recall /
%(fbeta^2*precision+recall); def = 1; the bigger the fbeta, the higher
%weight the recall is
%

if nargin < 3, correct_zerofeature = false; end
if nargin < 4, fbeta = 1; end

if nnz(size(answer)~=size(guess))> 0
    fprintf('answer ans guess are not of the same size\n');
    return
end

answer = full(answer);
guess = full(guess);
hit = sum(answer & guess, 1);
nguess = sum(guess, 1);
nanswer = sum(answer, 1);
precision = (hit ./ nguess)';
recall = (hit ./ nanswer)';

if correct_zerofeature
    precision(nanswer==0 & nguess==0) = 1;
    recall(nanswer==0 & nguess==0) = 1;
    precision(nanswer~=0 & nguess==0) = 0;
    recall(nanswer==0 & nguess~=0) = 0;
end

fmeasure = (1+fbeta^2)*(precision.*recall)./((fbeta^2)*precision+recall);
fmeasure(precision==0&recall==0) = 0;