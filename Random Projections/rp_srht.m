function P = rp_srht(n,s)
% Subsampled Hadamart transform. The creation method is not efficient!
% Auxiliary function for RP method in [1].
%
% This function computes the *normalized* subsampled hadamard transform
% matrix. 
%
% Inputs:
%   n   : the number of rows of the matrix
%   s   : the number of columns of the matrix
%
% Outputs:
%   P   : the random-projection matrix
%
% [1] E. Kontopoulou, G. Dexter, W. Szpankowski, A. Grama and P. Drineas, 
% "Randomized Linear Algebra Approaches to Estimate the von Neumann Entropy 
% of Density Matrices," in IEEE Transactions on Information Theory, 
% vol. 66, no. 8, pp. 5003-5021, Aug. 2020, doi: 10.1109/TIT.2020.2971991.
%
% Copyright: E. Kontopoulou, G. Dexter, W. Szpankowski, A. Grama, P. Drineas
%
% -- Last Update 11/06/2017


%-- Create S, the sampling matrix 
S= zeros(n,s);
j = randsample(n,s,true);
for i = 1:s
    S(j(i),i) = 1;
end

%-- Create H, the normalized Hadamard transform matrix
H = (1/sqrt(2)).*hadamard(n);

%-- Create D, the diagonal +-1 matrix
D = diag(randsample([-1,1],n,true));

%-- Construct P 
P = D*H*S;

