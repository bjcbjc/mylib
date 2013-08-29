function res = vbgm(data, K, varargin)
    %data: #samples x #features
    %
    %
    %

    [N, D]  = size(data);
    
    para.m_0 = zeros(D,1);
    para.b_0 = 0.005*N;
    para.a_0 = 0.01;
    para.nu_0 = D+1;
    para.W_0 = [];
    para.plotlb = false;
    para.maxloop = 1000;
    para.lbthres = 0.001;

    para = assignpara(para, varargin{:});
    if isempty(para.W_0)
        para.W_0 = eye(D);
%         para.W_0 = cov(data)\eye(D);
    end
           
    
    %initial responsibilities r
    fprintf('initial kmeans\n');
    assignments = kmeans(data, K);
    res.r = (ones(N,K) / K) * .1;
    for k = 1 : K
        res.r(:,k) = res.r(:,k) + .9 * (assignments == k);
    end

    fprintf('VB...\n');
    res.lb = NaN(para.maxloop,1);
    lower_bound = -Inf;
    do = 1;
    iter = 1;
    while do
        res = get_other_parameters(res, data, para, K, D);
        res = get_r(res, data, K, D, N);

        lower_bound_update = variational_lower_bound(res, data, para, K, D);        
        do = (lower_bound_update - lower_bound) > para.lbthres && iter <= para.maxloop ;
        if (lower_bound_update - lower_bound) < 0 && iter <= para.maxloop
            fprintf('decreasing lb, iter %d: %g -> %g\n',iter, lower_bound, lower_bound_update);
        end
        lower_bound = lower_bound_update;
        res.lb(iter) = lower_bound;
        iter = iter + 1;

        if para.plotlb
            plot(res.lb(1:iter-1));
            drawnow;
        else
            if mod(iter,100) == 0
                fprintf('%d iterations, lb =%g\n', iter, lower_bound);
            end
        end
    end
    res.lb(isnan(res.lb)) = [];
end

function paraest = get_other_parameters(paraest, X, para, K, D)
% This function calculated the updated parameter values for alpha, m, W, 
% nu, and beta given the value of r and the data matrix X.
%
%@param r           : n x k matrix for distribution of z_i
%@param X           : n x d data matrix
%
%@return alpha       : k x 1 matrix of positive dirichlet parameters
%@return m           : d x k matrix of means
%@return W           : k long cell array of d x d covariance matrics 
%@return nu          : k x 1 matrix of degrees of freedom for W matrices
%@return beta        : k x 1 matrix of scaling factors for NIW distributions

        
    Nk = sum(paraest.r)'; %K x 1
    xbar = bsxfun(@rdivide, paraest.r' * X, Nk); %k x D
    paraest.W = cell(K,1);
    for i = 1:K
        d = bsxfun(@minus, X, xbar(i,:)); %n x D        
        sk = d'*diag(paraest.r(:,i))*d; %without dividing Nk because later W{i} needs Nk*Sk_dividedNk
        d = xbar(i,:) - para.m_0'; %1 x D        
        paraest.W{i} = ( para.W_0\eye(D) + sk + (para.b_0*Nk(i)/(para.b_0+Nk(i))) * (d'*d) ) \ eye(D);
    end

    paraest.alpha = para.a_0 + Nk;
    paraest.beta = para.b_0 + Nk;
    paraest.nu = para.nu_0 + Nk + 1;
    paraest.m = bsxfun(@rdivide, ...
        bsxfun(@plus, para.b_0 * para.m_0', bsxfun(@times, Nk, xbar)), paraest.beta)'; %d x k

end

function paraest = get_r(paraest, X, K, D, N)
% This function calculates a matrix r which are parameters of the
% distribution over assignments for each row of the data matrix X.  The
% upated r is based on the values alpha, m, W, nu, and beta.
%
%@param alpha       : k x 1 matrix of positive dirichlet parameters
%@param m           : d x k matrix of means
%@param W           : k long cell array of d x d covariance matrics 
%@param nu          : k x 1 matrix of degrees of freedom for W matrices
%@param beta        : k x 1 matrix of scaling factors for NIW distributions
%@param X           : n x d data matrix
%
%@return r          : n x k matrix for distribution of z_i


    pik = psi(paraest.alpha) - psi(sum(paraest.alpha)) ;
    lambda = sum(psi(0.5*bsxfun(@minus, paraest.nu+1, 1:D)),2) + ...
        D*log(2) + log(max(realmin, cellfun(@det, paraest.W)));
    
    paraest.r = NaN(N, K);

    for j = 1:K
        d = bsxfun(@minus, X, paraest.m(:,j)');
        c1 = pik(j) + 0.5 * lambda(j);
        c2 = -0.5 * D/paraest.beta(j) ;
        for i = 1:N    
            paraest.r(i,j) = c1 + c2 - 0.5 * paraest.nu(j) * (d(i,:)*paraest.W{j}*d(i,:)') ;
        end
    end
    %-mean(paraest.r,2) to shift the values, avoiding Inf
    paraest.r = paraest.r - D/2*log(2*pi);    
%     if any(paraest.r > 500)
%         paraest.r = bsxfun(@minus, paraest.r, max(paraest.r,[],2)) + 500;
%     end
    %paraest.r = exp(paraest.r);
    paraest.r = exp( bsxfun(@minus, paraest.r, max(paraest.r, [], 2)) ); 
    %avoid r=0
%     paraest.r(paraest.r==0) = 1e-10;
%     paraest.r(isinf(paraest.r) &paraest.r>0) = realmax;
    paraest.r = bsxfun(@rdivide, paraest.r, sum(paraest.r,2));
%     paraest.r(paraest.r==0) = realmin;
end



function lb = variational_lower_bound(paraest, X, para, K, D)
% This function calculates the variational lower bound.  This value should
% go up as the algorithm progresses.  The algorithm only stops when the
% variational lower bound has stopped increasing.
%
%@param r           : n x k matrix for distribution of z_i
%@param alpha       : k x 1 matrix of positive dirichlet parameters
%@param m           : d x k matrix of means
%@param W           : k long cell array of d x d covariance matrics 
%@param nu          : k x 1 matrix of degrees of freedom for W matrices
%@param beta        : k x 1 matrix of scaling factors for NIW distributions
%@param X           : n x d data matrix
%
%@return lb         : calculated scalar lower bound


    logpik = psi(paraest.alpha) - psi(sum(paraest.alpha)); %k x 1
    loglambdak = D*log(2) + sum(psi(0.5 * bsxfun(@minus, paraest.nu+1, (1:D))),2) + ...
        log(max(realmin, cellfun(@det, paraest.W))); %k x 1
    Nk = sum(paraest.r)';
    xbar = bsxfun(@rdivide, paraest.r' * X, Nk); % k x D
    Sk = cell(K, 1);
    ePX_ZuL = 0;
    ePuL = 0;
    eQuL = 0;
    for i = 1:K
        d = bsxfun(@minus, X, xbar(i,:)); %n x D
        Sk{i} = d'*diag(paraest.r(:,i))*d ./ Nk(i);
        dk = xbar(i,:) - paraest.m(:,i)'; %1 x D
        ePX_ZuL = ePX_ZuL + Nk(i) * ( loglambdak(i) - D/paraest.beta(i) ...
            - paraest.nu(i) * trace( Sk{i} * paraest.W{i} ) ...
            - paraest.nu(i) * (dk*paraest.W{i}*dk') - D*log(2*pi) );
        dk = (paraest.m(:,i) - para.m_0)'; %1 x D
        ePuL = ePuL + D * log(para.b_0/2/pi) + loglambdak(i) - ...
            D*para.b_0/paraest.beta(i) - para.b_0*paraest.nu(i)* (dk*paraest.W{i}*dk');        
        eQuL = eQuL + (paraest.nu(i)-D)*loglambdak(i) + ...
            D*(log(paraest.beta(i)/2/pi) -1 -paraest.nu(i)) + 2*lnB(paraest.W{i}, paraest.nu(i), D);
    end

    invW0 = para.W_0 \ eye(D);
    ePX_ZuL = 0.5 * ePX_ZuL;
    ePuL = 0.5*ePuL + K*lnB(para.W_0, para.nu_0, D) + 0.5*(para.nu_0-D-1)*sum(loglambdak) ...
        - 0.5*sum( paraest.nu .* cellfun(@(wk) trace(invW0*wk), paraest.W) );
    eQuL = 0.5*eQuL;

    ePZ_pi = sum(sum( bsxfun(@times, paraest.r, logpik') ) );
    ePpi = gammaln(para.a_0*K) - K*gammaln(para.a_0) + (para.a_0-1) * sum(logpik);
    eQZ = sum(sum( paraest.r .* log(paraest.r) ));
    eQpi = sum((paraest.alpha-1) .* logpik) + gammaln(sum(paraest.alpha)) - sum(gammaln(paraest.alpha));
    
    lb = ePX_ZuL + ePZ_pi + ePpi + ePuL - eQZ - eQpi - eQuL;

end

function lB = lnB(W, nu, D)
    lB = -0.5*nu* (log(max(realmin,det(W))) + D*log(2)) ...
        - 0.25*D*(D-1)*log(pi) - sum(gammaln( (nu+1-(1:D)')/2 ));
end
