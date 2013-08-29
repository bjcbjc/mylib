function [K, H, G] = canonGaussian(mu, invCov)
n = length(mu);
H = mu*invCov;
K = invCov;
% G = -.5*mu'*invCov*mu- 0.5*log( det(invCov\(eye(n))) ) - (n/2) * log(2*pi)   ;
G = [];