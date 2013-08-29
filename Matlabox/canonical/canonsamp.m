function sampX = canonsamp(distr)
    %sample n data from distr
    %
    %distr: canonical form of the distribution, multivariate gaussian
    %n: number of data points the user wishes to sample
    %
    %sampX: sample data, d x n, d is the number of variables, ordered the
    %   same as distr.scope
    %   

    d = length(distr.scope);
    
    %mvnrnd: sigma is covariance
    mu = distr.k \ distr.h;
    sigma = distr.k \ eye(d);
    
    sigma = errcheck(sigma);
    
    sampX = mvnrnd(mu', sigma); %n x d
    
end

function sigma = errcheck(sigma)
    [~, err] = cholcov(sigma);
    if err ~= 0 %numeric precision
        if ~all(all(sigma==sigma'))
            sigma = max(sigma, sigma');
        else
            [W V] = eig(sigma);
            d = diag(V);
            if any(d<0)
                d(d<0) = 0;
                V = diag(d);
                sigma = W * (V) * W';
            else
                sigma = W * (V+eps) * W';
            end
        end
    end
    [~, err] = cholcov(sigma);
    if err ~= 0 %numeric precision        
        [W V] = eig(sigma);
        d = diag(V);
        if any(d<0)
            d(d<0) = 0;
            V = diag(d);
            sigma = W * (V) * W';
        else
            sigma = W * (V+eps) * W';
        end        
    end
end