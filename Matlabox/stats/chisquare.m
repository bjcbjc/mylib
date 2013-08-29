function [s] = chisquare(table)
%calculate chisq statistics
%input:
%table: 2x3 or 2x2
%return:
%s: chi-square statistices
%
    colsum = sum(table);
    rowsum = sum(table')';    
    n = sum(colsum);
    pj = colsum ./ n;
    pi = rowsum ./ n;
    df = numel(table)-1;
    
    E = n * (pi * pj);
    diff = abs(table - E);
    
    %Yate's correction if any count is less than 5
    [i j] = find(table<5);
    if ~isempty(i) && numel(table)==4
        diff = diff - 0.5;
        %fprintf('yates\n');
    end
%     for index = 1:length(i)
%         diff(i(index),j(index)) = diff(i(index),j(index)) - 0.5;
%     end
    %s = nansum(nansum(((table - E) .^ 2) ./ E));
    s = nansum(nansum((diff .^ 2) ./ E));
end