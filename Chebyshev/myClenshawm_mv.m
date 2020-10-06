function y = myClenshawm_mv(R,m,a,u,g)
% This Clenshaw implementation supports multiplication with
% gaussian vectors. Auxiliary function for Chebyshev method appeared in [1]
%
% Inputs:
%
% R       : n x n input matrix (density matrix)
% m       : number of Chebyshev summands to keep (counts from 0)
% a       : Chebyshev coefficients (this should be m+1)
% g       : Gaussian vector
% 
% Ouput:
%
% Y       : the value of the function on R.
%
% [1] E. Kontopoulou, G. Dexter, W. Szpankowski, A. Grama and P. Drineas, 
% "Randomized Linear Algebra Approaches to Estimate the von Neumann Entropy 
% of Density Matrices," in IEEE Transactions on Information Theory, 
% vol. 66, no. 8, pp. 5003-5021, Aug. 2020, doi: 10.1109/TIT.2020.2971991.
%
% Copyright: E. Kontopoulou, G. Dexter, W. Szpankowski, A. Grama, P. Drineas
%
% Last Update: 10/28/2017

n = size(R,2);

y1 =zeros(n,1);
y2 = zeros(n,1);
y0 = zeros(n,1);

for k = m+1:-1:1
    y2 = y1;
    y1 = y0;
    y0 = (4/u).*(R*y1) - (2.*y1 + y2) + a(k).*g ;
end
y = (1/2)*((g'*g)*a(1) + g'*(y0 - y2));