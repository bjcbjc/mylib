function [res] = mic_single_mtx(X, Y, W, ntoptest)
    % Initialize values.

    if nargin < 4
        ntoptest = 100;
    end
    
    if nargin < 3
        W = ones(size(Y));
    end
    
    [n,p] = size(X);
    [n2,nY] = size(Y);    
    assert(n==n2, 'number of observation is not the same');

    costPerCoeff = 2;
    pForFeatureCost = p;

    % Create a 'model' struct to store important state.
    model.X = X;
    model.Y = Y;
    model.W = sqrt(W);

    % Store a method that, given residuals, computes sigmaHat.
    model.sigmaHatOfResiduals = @(x) var(x); %%%%% implement to be diag var (residual)

    % for each added feature
    model.costPerFeature = log2(pForFeatureCost) + costPerCoeff;

    % Store info on what features we've used.
    model.featurepool = 2:p;
    
    global nrmX
    if all(model.W(:)==1)
        nrmX = zeromean_univar_normalization(X,1);
        model.weighted = false;
    else
        model.weighted = true;
    end
    
    model.betaHat = sparse(p,nY);
    
    for yi = 1:nY
        %%%%%%%%%%%%%%%%%%%%%%%%%% reset these attributes %%%%%%%%%%%%%%%%%%
        %initialize this for every yi, because a marker is added, the cost
        %is set to 0
        % available features
        model.featurepool = 2:p;

        % Start computing betaHat.        
        model.betaHat(1,yi) = mean(model.Y(:,yi)); % Initial intercept terms.
        model.curResiduals = model.Y(:,yi) - model.X * model.betaHat(:,yi);
        model.curSigmaHat = model.sigmaHatOfResiduals(model.curResiduals);
        
        model.likelihoodCost = likeCostAdd(1,model,yi);
        model.codingCost = 0; % no cost in the beginning, even though intercept is used
        
        bestCostSoFar = model.likelihoodCost;
        if model.weighted
            nrmX = zeromean_univar_normalization(bsxfun(@times, model.W(:,yi), X), 1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Try adding features until it no longer helps.
        for nFeatures = 2:p
            %take the feature that increases the likelihood most and see of the
            %cost is worthy
            %[likeCosts candBetaHats] = arrayfun(@(f)likeCostOfAddFeature(f,model), ...
            %    model.featurepool, 'UniformOutput',false);
            %likeCosts = cell2mat(likeCosts);
            %[tmp toAdd] = min(likeCosts);

            [likeCost candBetaHat toAdd] = bestlikeCostAdd(model, yi, model.featurepool, ntoptest);

            costIfAdded = likeCost + model.codingCost + model.costPerFeature;

            if costIfAdded < bestCostSoFar
%                 bestCostSoFar = costIfAdded;

                % Add that feature to the model.
                % Update state.
                model.betaHat(:,yi) = candBetaHat;
                model.likelihoodCost = likeCost;
                model.codingCost = model.codingCost + model.costPerFeature;
                model.curResiduals = model.Y(:,yi) - model.X * model.betaHat(:,yi);
                model.curSigmaHat = model.sigmaHatOfResiduals(model.curResiduals);
                model.likeforcmp = negloglike(model, model.curResiduals, yi);
                bestCostSoFar = model.likeforcmp + model.codingCost;
                model.featurepool(toAdd) = [];
            else
                break;
            end
        end
    end

    % Return betaHat.
    res.beta = model.betaHat;

    % Prepare return results.
    res.yhat = X * model.betaHat;
    res.weighted = model.weighted;
end

function [negLogLike candidateNewBetaHat] = likeCostAdd(feature,model,yi)

    oldBetaHat = model.betaHat(:,yi);
    candidateNewBetaHat = oldBetaHat;   
    if(oldBetaHat(feature)==0)
        allf = union(feature,find(oldBetaHat));
        candidateNewBetaHat(allf) = ...
            bsxfun(@times, model.W(:,yi), model.X(:,allf)) \ ...
            bsxfun(@times, model.W(:,yi), model.Y(:,yi));
    end
    
    % Get residuals.
    epsHat = model.Y(:,yi) - model.X * candidateNewBetaHat;
    
    %why the default is not to update sigmaHat?
    %because when calcualte the likelihood, sigma should be the previoud
    %one, ie. the one before the new feature is added
    %sigmaHat = model.curSigmaHat;

    % negLogLikelihood
%     [n h] = size(epsHat);
%     firstTerm = (n*h/2) * log2(2*pi) + (n/2)*log2(model.curSigmaHat);
%     undividedSecondTerm = sum(epsHat.^2) / model.curSigmaHat / (2*log(2));
%     negLogLike = firstTerm + undividedSecondTerm;
    negLogLike = negloglike(model, epsHat, yi);
end

function [negLogLike newBetaHat maxi] = bestlikeCostAdd(model, yi, testfi, ntoptest)
    global nrmX 

    %prepare selected X
    tmpX = bsxfun(@times, model.W(:,yi), model.X(:,model.betaHat(:,yi)~=0));

    %calculate corr
    nrm_residual = zeromean_univar_normalization(...
        bsxfun(@times, model.W(:,yi), model.curResiduals), 1);
    
    c = abs( singOlsMtx(nrmX(:,testfi), nrm_residual));
    
    [~, si] = sort(c, 'descend');

    wY = bsxfun(@times, model.W(:,yi), model.Y(:,yi));
    wsigma = model.curSigmaHat ./ (model.W(:,yi).^2);
    %take the top corr to calculate the likelihood
    tmpeps = ones(min([ntoptest length(testfi)]),1)*Inf;
    for i = 1:min([ntoptest length(testfi)])
        toAddi = testfi(si(i));
        toAddX = bsxfun(@times, model.W(:,yi), model.X(:,toAddi));
        tmpeps(i) = sum( ((wY - ...
            [tmpX toAddX] * ([tmpX toAddX] \ wY)).^2) ./ wsigma);
    end
        
    [~, maxi] = min(tmpeps);
    maxi = si(maxi);    

    toAddi = testfi(maxi);
    toAddX = bsxfun(@times, model.W(:,yi), model.X(:,toAddi));
        
    b = [tmpX toAddX] \ wY;
    
    %update betaHat
    newBetaHat = model.betaHat(:, yi); %don't take the whole matrix, save some memory    
    newBetaHat(model.betaHat(:,yi)~=0) = b(1:end-1);    
    newBetaHat( toAddi ) = b(end);

    epsHat = model.Y(:,yi) - [tmpX toAddX] * b; %only for yi
    
    % negLogLikelihood
%     [n nY] = size(epsHat);
%     firstTerm = (n*nY/2) * log2(2*pi) + (n/2)*log2(model.curSigmaHat);
%     secondTerm = sum(epsHat.^2) / model.curSigmaHat / (2*log(2));
%     negLogLike = firstTerm + secondTerm;
    negLogLike = negloglike(model, epsHat, yi);
end

function nloglike = negloglike(model, epsHat, yi)
    n = size(epsHat,1);
    wsigma = model.curSigmaHat ./ (model.W(:,yi).^2);
    firstTerm = (n/2) * log2(2*pi) + (1/2)*sum(log2(wsigma));
    undividedSecondTerm = sum((epsHat.^2)./wsigma) / (2*log(2));
    nloglike = firstTerm + undividedSecondTerm;
end