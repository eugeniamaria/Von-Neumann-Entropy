function H = ChebEntropy(R,m,p,u)
% This function computes the Entropy of the density matrix R, using the
% Chebyshev approximation of the first kind and the **Gaussian** trace 
% estimator. Method appeared in [1].
%
% Inputs:
%
% R       : n x n input matrix (density matrix)
% m       : number of Chebyshev summands to keep
% p       : number of Gaussian vectors
% u       : the estimation of the largest eigenvalue
% 
% Ouput:
%
% H       : the VonNeumann Entropy 
%
% [1] E. Kontopoulou, G. Dexter, W. Szpankowski, A. Grama and P. Drineas, 
% "Randomized Linear Algebra Approaches to Estimate the von Neumann Entropy 
% of Density Matrices," in IEEE Transactions on Information Theory, 
% vol. 66, no. 8, pp. 5003-5021, Aug. 2020, doi: 10.1109/TIT.2020.2971991.
%
% Copyright: E. Kontopoulou, G. Dexter, W. Szpankowski, A. Grama, P. Drineas
%
% Last Update: 10/29/2017

n = size(R,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u = min([6*(powermethod(R,ceil(log(n)),ceil(log(1/delta)))),1]); 
% estimate the largest eigenvalue
% u = min([6*eigs(R,1);1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = randn(n,p); % Gaussian vectors
y = zeros(p,1);
%--Compute the chebyshev series
%   -- step 1: construct the coefficients
a = chebyshevcoeff(m,u);
%   -- step 2: compute the chebyshev polynomial on R
for j=1:p
   y(j) = myClenshawm_mv(R,m,a,u,G(:,j));
end

H = -(1/p) * sum(y);

fprintf('\t Entropy = %f\n', H);