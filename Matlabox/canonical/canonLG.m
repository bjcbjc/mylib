function [K, H, G] = canonLG(W, beta)
%This function expresses the linear gaussian conditional factor 
%characterized by W and beta in canonical form according to the equations
%specified in the technical report: Implementation of Clique Tree Algorithm 
%for Conditional Linear Gaussian Model by Wei Sun.

%W is the coefficient vector and beta is the precision scalar.  W must have
%at least two elements since it is assumed to contain the constant
%coefficient as the first index.  W is supposed to be a column vector.  In 
%the output, the first index of H (and the first row/column of K correspond
%to the child variable in the regression model while subsequent indices 
%occur in the order of the parent variables corresponding to W.
mu = W(1);
W = W(2:end);

H = [beta*mu;-beta*mu.*W];
K = [beta, -beta.*W';-beta.*W,beta.*W*W'];
% G = -.5*((mu^2)*beta - log(2*pi*beta));
G = [];

