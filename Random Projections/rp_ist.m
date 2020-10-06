function P = rp_ist(n,s)
% This function computes the Input-Sparcity Transform matrix.
% Auxiliary function for RP method in [1].
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
j = randsample(s,n,true);
for i = 1:n
    S(i,j(i)) = 1;
end

%-- Create D, the diagonal +-1 matrix
D = diag(randsample([-1,1],n,true));

%-- Construct P 
P = D*S;