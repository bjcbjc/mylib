classdef mulViewCluster < handle
    
    properties
        X;
        
        %hyper-para
        alpha;
        beta;
        mean0;
        cov0;
        psi0;
        p0;
        
        %fake infinity
        nCluster;
        nView;
        
        %cluster parameteres
        Xbar;
        Mean;
        Cov;        
        
        %mixture; posterior
        Pi;
        Eta;        
        
        %parameteres for CRP
        betaParaForPi;
        betaParaForEta;
        
        %latent variable: assignment
        Y;
        Z;
    end
    
    
    methods
        function obj = mulViewCluster(X, varargin)
            obj.X = X;
            obj.alpha = 1;
            obj.beta = 1;            
            [N, F] = size(X);
            obj.mean0 = [];
            obj.cov0 = [];
            obj.psi0 = [];
            obj.p0 = [];
            obj.nCluster = N;
            obj.nView = F;
            
            obj = assignpara(obj, varargin{:});
            obj.nCluster = min(N, obj.nCluster);
            obj.nView = min(F, obj.nView);
            
            if isempty(obj.mean0)
                obj.mean0 = zeros(1, F);
            end
            if isempty(obj.cov0)
                obj.cov0 = eye(F);
            end
            if isempty(obj.psi0)
                obj.psi0 = 1;%ones(obj.nView, obj.nCluster);
            end
            if isempty(obj.p0)
                obj.p0 = 1; %ones(obj.nView, obj.nCluster);
            end
            
            obj.Xbar = cell(obj.nView, 1);
            obj.Sigma = cell(obj.nView, obj.nCluster);
            obj.Mean = cell(obj.nView, 1);
            obj.Cov = cell(obj.nView, obj.nCluster);
            for v = 1:obj.nView
                obj.Xbar{v} = zeros(obj.nCluster, F);
                obj.Mean{v} = zeros(obj.nCluster, F);
                for k = 1:obj.nCluster
                    obj.Cov{v,k} = eye(F);
                    obj.Sigma{v,k} = eye(F);
                end
            end
            
            obj.Pi = ones(F, obj.nView);
            obj.Eta = ones(obj.nView, N, obj.nCluster);
            obj.Y = zeros(F, obj.nView); %binary
            obj.Z = zeros(obj.nView, N); %multinomial 1-obj.nCluster           
        end
        
        function initialize(obj)
        end
        
        function updateDistribution(obj)
            
            for v = 1:obj.nView
                for k = 1:obj.nCluster
                    obj.Xbar{v}(k, :) = 0;
                    validFeature = obj.Y(:, v) == 1;
                    members = obj.Z(v,:) == k;
                    obj.Xbar{v}(k, validFeature) = nanmean( obj.X( members, validFeature), 1);
                    
                    obj.Mean{v}(k, :) = 0;
                    n_vk = nnz( members );
                    obj.Mean{v}(k, validFeature) = ( n_vk * obj.Xbar{v}(k, validFeature) + ...
                        obj.psi0 * obj.mean0(validFeature) ) ./ (n_vk + obj.psi0);
                    
                    obj.Cov{v,k}(:) = 0;
                    distToMean = bsxfun(@minus, obj.X( members, validFeature), obj.Xbar{v}(k, validFeature));
                    distToMean0 = bsxfun(@minus, obj.X( members, validFeature), obj.mean0(validFeature));
                    obj.Cov{v,k}(validFeature, validFeature) = obj.cov0(validFeature, validFeature) + ...
                        distToMean' * distToMean + ...
                        n_vk * obj.psi0 * (distToMean0' * distToMean0) ./ (n_vk + obj.psi0);                    
                end
            end
                
        end
        
        function updateAssignment(obj)
            PiPortion = zeros(F, obj.nView);
            obj.betaParaForPi(:, 1) = 1 + sum(obj.Pi,1); %1 x nView            
            for v = 1:obj.nView-1                
                obj.betaParaForPi(v, 2) = obj.alpha + sum(sum(obj.Pi(:, v+1:end)));            
                obj.betaParaForEta(v, :, 1) = 1 + sum(squeeze(obj.Eta(v, :, :)), 2);
                for k = 1:obj.nCluster-1
                    obj.betaParaForEta(v, k, 2) = obj.beta + sum(sum(obj.Eta(v, :, k+1:end)));
                end
                logPostMultPi = psi(obj.betaParaForPi(v,1)) - psi(obj.betaParaForPi(v,1)+obj.betaParaForPi(v,2)) + ...
                    sum( psi( obj.betaParaForPi(1:v-1,2)) - psi( sum(obj.betaParaForPi(1:v-1,:),2) )  );
                lq = obj.logQ(v);                
                PiPortion(obj.Y(:, v) == 1, v) = logPostMultPi + lq;
            end
            obj.Pi
        end
        
        function lq = logQ(obj, v)
            validFeature = obj.Y(:, v) == 1;
            n_v = obj.getNumFeatureInClusterPerView(v);
            d_v = sum( obj.Y(:, v) == 1);
            d_v2 = d_v/2;
            psi_v = obj.psi0 + n_v; %vector of K
            p_v = obj.p0 + n_v; %vector of K
            lq = 0;
            for k = 1:obj.nCluster
                memberIdx = find(obj.Z(v, :) == k);
                invSigma = pinv(obj.Sigma{v,k}(validFeature, validFeature) );
                logDetSigma = sum( psi((p_v(k)+1-(1:d_v(k)))/2) ) + d_v(k)*log(2) + log( det( obj.Cov{v,k}(validFeature, validFeature)));  
                lq = lq + n_v(k) * ( ...
                    d_v2 * log( psi_v(k) ) - (dv+p_v(k)) * log(2) - (d_v*d_v/4 + d_v/4)*log(pi) + ...
                    (p_v(k)/2 + d_v2 + 1) * logDetSigma + ...
                    p_v(k)/2 * log( det( obj.Cov{v,k}(validFeature, validFeature) ) ) - ...
                    0.5 * trace( obj.Cov{v,k}(validFeature, validFeature) * invSigma ) - ...
                    sum( gammaln( (p_v(k)+1-(1:dv))/2 ) ) );
                
                for xi = memberIdx'
                    lq = lq - 0.5 * obj.X(xi, validFeature) * (psi_v(k) * invSigma) * obj.X(xi, validFeature)';
                end                    
            end            
        end
        
        function n_v = getNumFeatureInClusterPerView(obj, v)
            n_v = zeros( obj.nCluster, 1);
            for k = 1:obj.nCluster
                n_v(k) = sum( obj.Z(v, :) == k );
            end
        end
    end
    
    
end


