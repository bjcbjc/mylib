function gseaplot(score, annotation, func, P, fontsize)
    if nargin < 3,
        func = @abs;
    end
    if nargin < 4,
        P = 1;
    end
    if nargin < 5
        fontsize = 12;
    end


    if size(score,2) > size(score,1)
        score = score';
    end
    
    if size(annotation,2) > size(annotation,1)
        annotation = annotation';
    end


    if size(annotation,2) ~= 1 || size(score,2) ~= 1
        fprintf('gseaplot only does one annoation at a time\n');
        return
    end

    penalty = -1/sum(~annotation);
    transformed_score = func(score);
    enrichscore = transformed_score .^ P;
    NR = sum(bsxfun(@times, abs(transformed_score), annotation) .^ P, 1);
    enrichscore = bsxfun(@rdivide, enrichscore, NR);
    enrichscore(~annotation) = penalty;
    
    [~, si] = sort(score, 'descend');
    
    n = size(score,1);
    tickunit = round(n/5/1000)*1000;
    
    subplot(5,1,1)
    nanimagesc( annotation(si)', genColorMap('kw', 2))
    set(gca, 'ytick', [], 'xtick', 0:tickunit:n, 'fontsize', fontsize-2)
    xlim([0, n+0.5])
    title('annotation  ', 'fontsize', fontsize)
    
    subplot(5,1,2)
    h = area( score(si)');
    set(gca, 'xtick', 0:tickunit:n, 'fontsize', fontsize-2)
    set(h, 'facecolor', [0.1 0.7 0.5], 'linestyle', '-');
    ylabel('score  ', 'fontsize', fontsize)
    xlim([0, n+0.5])
    
    subplot(5,1,3:5)
    plot(cumsum(enrichscore(si)), 'b.-');
    set(gca, 'xtick', 0:tickunit:n, 'fontsize', fontsize-2)
    ylabel('ES  ', 'fontsize', fontsize)
    xlim([0, n+0.5])
    
    