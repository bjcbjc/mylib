function [res] = micregress(X, Y, varargin)
    %
    %X: #sample x #feature, the first column has to be all one (intercept)
    %Y: #sample x 1
    %
    %options:
    %   ntoptest: number of top correlated features for testing (for speed
    %       up); set it to 0 if no speed-up (all features will be tested)
    %   costPerCoeff: cost of each coefficient, default 2
    %   costprior: default, uniform prior which specifies the cost of each
    %       feature = log2(#features); if prior is specified, it should be
    %       a vector of #feature element, in which the cost is
    %       -log2(probability of the feature should be chosen)
    %
    
    para.ntoptest = 100;
    para.costPerCoeff = 2;
    para.costprior = log2(size(X,2));
    para.sigmathres = 1e-8;
    para.ridge = 0;
    
    para = assignpara(para, varargin{:});
    
    [n,p] = size(X);
    [n2,nY] = size(Y);    
    assert(n==n2, 'number of observation is not the same');
    assert(nY==1, 'only support one Y now');
    
    % Create a 'model' struct to store important state.
    model.X = X;
    model.Y = Y;
    model.n = n;
    
    % Store a method that, given residuals, computes sigmaHat.
    model.sigmaHatOfResiduals = @(x) var(x); %%%%% implement to be diag var (residual)

    % for each added feature
    model.costPerFeature = para.costprior + para.costPerCoeff;
    if length(model.costPerFeature) == 1 %uniform
        model.costPerFeature = repmat(model.costPerFeature, p, 1);
    end
    if size(model.costPerFeature, 2) > size(model.costPerFeature, 1)
        model.costPerFeature = model.costPerFeature';
    end

    % Store info on what features we've used.
    model.featurepool = 2:p;
    
    
    global nrmX
    
    if para.ntoptest > 0
        nrmX = zeromean_univar_normalization(X,1);
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
        
        model.likelihoodCost = likeCostAdd(1,model,yi, para);
        model.codingCost = 0; % no cost in the beginning, even though intercept is used
        
        bestCostSoFar = model.likelihoodCost;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Try adding features until it no longer helps.
        for nFeatures = 2:p
            %take the feature that increases the likelihood most and see of the
            %cost is worthy
            
            [likeCost candBetaHat toAdd addfi] = bestCostAdd(model, yi, model.featurepool, para);

            costIfAdded = likeCost + model.codingCost + model.costPerFeature(addfi);

            if costIfAdded < bestCostSoFar
%                 bestCostSoFar = costIfAdded;

                % Add that feature to the model.
                % Update state.
                model.betaHat(:,yi) = candBetaHat;
                model.likelihoodCost = likeCost;
                model.codingCost = model.codingCost + model.costPerFeature(addfi);
                model.curResiduals = model.Y(:,yi) - model.X * model.betaHat(:,yi);
                model.curSigmaHat = model.sigmaHatOfResiduals(model.curResiduals);
                model.likecostforcmp = negloglike(model.n, model.curSigmaHat, sum(model.curResiduals.^2));
                model.featurepool(toAdd) = [];
                bestCostSoFar = model.likecostforcmp + model.codingCost;
                if model.curSigmaHat < para.sigmathres
                    break;
                end
            else
                break;
            end
        end
    end

    % Return betaHat.
    res.beta = model.betaHat;

    % Prepare return results.
    res.yhat = X * model.betaHat;
    res.sigmahat = model.curSigmaHat;
end

function [negLogLike candidateNewBetaHat] = likeCostAdd(feature,model,yi, para)

    oldBetaHat = model.betaHat(:,yi);
    candidateNewBetaHat = oldBetaHat;   
    if(oldBetaHat(feature)==0)
        allf = union(feature,find(oldBetaHat));
        if para.ridge ~= 0
            candidateNewBetaHat(allf) = semiridgeregress(model.X(:,allf), model.Y(:,yi), para.ridge);
        else
            candidateNewBetaHat(allf) = model.X(:,allf) \ model.Y(:,yi);
        end
    end
    
    % Get residuals.
    epsHat = model.Y(:,yi) - model.X * candidateNewBetaHat;
    
    %why the default is not to update sigmaHat?
    %because when calcualte the likelihood, sigma should be the previoud
    %one, ie. the one before the new feature is added
    %sigmaHat = model.curSigmaHat;

    % negLogLikelihood
    negLogLike = negloglike(model.n, model.curSigmaHat, sum(epsHat.^2));
end

function [negLogLike newBetaHat maxi toAddi] = bestCostAdd(model, yi, testfi, para)
    global nrmX 
    
    
    %prepare selected X
    tmpX = model.X(:, model.betaHat(:,yi)~=0);

    if para.ntoptest > 0 %only screen for ntoptest correlated features
        %calculate corr
        nrm_residual = zeromean_univar_normalization(model.curResiduals, 1);
        
        c = abs( fastcorr(nrmX(:,testfi), nrm_residual, model.n) );
        
        [~, si] = sort(c, 'descend');
                
        %take the top corr to calculate the likelihood
        ntest = min([para.ntoptest length(testfi)]);
        tmpeps = ones(ntest, 1) * Inf;
        for i = 1:ntest
            toAddi = testfi(si(i));
            toAddX = model.X(:,toAddi);
            %sum squares of residuals
            if para.ridge ~= 0
                tmpb = semiridgeregress([tmpX toAddX], model.Y(:,yi), para.ridge);                
                tmpeps(i) = sum( ( model.Y(:,yi) - ...
                    [tmpX toAddX] * tmpb ).^2  );
            else
                tmpeps(i) = sum( ( model.Y(:,yi) - ...
                    [tmpX toAddX] * ([tmpX toAddX] \ model.Y(:,yi)) ).^2  );
            end
        end
        
        testcost = negloglike(model.n, model.curSigmaHat, tmpeps);
        %the other part of coding cost is the same; so only need to compare
        %the new likecost + the cost of the new feature
        testcost = testcost + model.costPerFeature( testfi( si(1:ntest) ) );
        
        [~, maxi] = min(testcost);
        maxi = si(maxi);
               
        
    else %screen for all available features
        %sum squares of residuals
        tmpeps = ones(length(testfi), 1)*Inf;
        for i = 1:length(testfi)
            toAddi = testfi(i);
            toAddX = model.X(:,toAddi);
            if para.ridge ~= 0
                tmpb = semiridgeregress([tmpX toAddX], model.Y(:,yi), para.ridge);
                tmpeps(i) = sum( ( model.Y(:, yi) - ...
                    [tmpX toAddX] * tmpb ).^2 );
            else
                tmpeps(i) = sum( ( model.Y(:, yi) - ...
                    [tmpX toAddX] * ([tmpX toAddX] \ model.Y(:,yi)) ).^2 );
            end
        end
        
        testcost = negloglike(model.n, model.curSigmaHat, tmpeps);
        %the other part of coding cost is the same; so only need to compare
        %the new likecost + the cost of the new feature
        testcost = testcost + model.costPerFeature( testfi );
        
        [~, maxi] = min(testcost);
    end        
    

    toAddi = testfi(maxi);
    toAddX = model.X(:,toAddi);
        
    if para.ridge ~= 0
        b = semiridgeregress([tmpX toAddX], model.Y(:, yi), para.ridge);
    else
        b = [tmpX toAddX] \ model.Y(:, yi);
    end
    
    %update betaHat
    newBetaHat = model.betaHat(:, yi); %don't take the whole matrix, save some memory    
    newBetaHat(model.betaHat(:,yi)~=0) = b(1:end-1);    
    newBetaHat( toAddi ) = b(end);

    epsHat = model.Y(:,yi) - [tmpX toAddX] * b; %only for yi
    
    % negLogLikelihood
    negLogLike = negloglike(model.n, model.curSigmaHat, sum(epsHat.^2));
end


function like = negloglike(n, sigma, ssqres)    
    
    %residual can be a matrix, each row is a tested model and each column
    %is a data point
    %negloglike returns -loglike for size(residual,1) models
                
    nmodel = size(ssqres,1);    
    if nmodel > 1 && size(sigma, 1) == 1
        sigma = repmat(sigma, nmodel, 1);
    elseif nmodel > 1 && nmodel ~= size(sigma, 1)
        error('negloglike: nmodel inconst in sigma and residual\n');
    end
    
%     ssqres = sum( residual.^2, 2);
    
    vec1 = (n/2) * log2(2*pi) + (n/2).*log2(sigma); % #nmodel x 1
    vec2 = ssqres./sigma ./ (2*log(2)); % #nmodel x #1
    like = vec1 + vec2; % #test-model x 1     
end

