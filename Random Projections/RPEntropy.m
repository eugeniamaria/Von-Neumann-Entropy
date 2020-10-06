function [H,time] = RPEntropy(R,k,s,method)
% Random Projection Entropy 
%
% This function computes the entropy of a rank deficient Density matrix
% using random projections. The random projection matrices can be
% constructed using two different methods:
%  [1] The normalized subsampled Hadamard Transform 
%  [2] The Input-Sparsity Transform
%  [3] Gaussian Random Projection
% Method appeared in [1].
%
% Inputs:
%
% R       : the density matrix
% k       : the target rank
% s       : the columns of the random projection matrix
% method  : the type of projection matrix : 'sRHT' for subsampled randomized
%           Hadmard Transform, 'IST' for Input-Sparsity Transform, 'GP' for
%           Gaussian projection.
% 
% Ouput:
%
% H      : the entropy of the matrix
%
% [1] E. Kontopoulou, G. Dexter, W. Szpankowski, A. Grama and P. Drineas, 
% "Randomized Linear Algebra Approaches to Estimate the von Neumann Entropy 
% of Density Matrices," in IEEE Transactions on Information Theory, 
% vol. 66, no. 8, pp. 5003-5021, Aug. 2020, doi: 10.1109/TIT.2020.2971991.
%
% Copyright: E. Kontopoulou, G. Dexter, W. Szpankowski, A. Grama, P. Drineas
%
% -- Last Update 10/06/2020

n = size(R,2);

disp('Constructing Random Projection')
%--Step 1: Construction of the Random Projection Matrix
if strcmp(method,'sRHT')
    %-- Subsampled Randomized Hadamard Transform (sRHT)
    tic,P = rp_srht(n,s); contime = toc;
    fprintf('\t Random Projection construction = %f sec\n', contime);
elseif strcmp(method,'IST')
    %-- Input-Sparsity Transform (IST)
    tic,P = rp_ist(n,s); contime = toc;
    fprintf('\t Random Projection construction = %f sec\n', contime);
elseif strcmp(method,'GP')
    %-- Gaussian Projection(GP)
    tic,P = randn(n,s); contime = toc;
    fprintf('\t Random Projection construction = %f sec\n', contime);
else
    error('Unknown projection method.')
end

tic,
%--Step 2: Projection
Rhat = R*P;

%--Step 3: Computation of the entropy
l = svds(Rhat,k);
sv = l./sum(l);
H = RealEntropy(sv);
time = toc;