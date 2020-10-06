function a = chebyshevcoeff(m,u)
% This function computes the coefficients of the shifted chebyshev
% polynomial of the first kind theat approximates the function f(x) =
% xlogx. Auxiliary function for Chebyshev method appeared in [1].
%
% Inputs:
%
%   m : degree of the polynomial 
%   u : upper bound of the interval that the polynomial exists
%
% Output:
%
%   a : the coefficients of the shifted chebyshev polynomial
%
% [1] E. Kontopoulou, G. Dexter, W. Szpankowski, A. Grama and P. Drineas, 
% "Randomized Linear Algebra Approaches to Estimate the von Neumann Entropy 
% of Density Matrices," in IEEE Transactions on Information Theory, 
% vol. 66, no. 8, pp. 5003-5021, Aug. 2020, doi: 10.1109/TIT.2020.2971991.
%
% Copyright: E. Kontopoulou, G. Dexter, W. Szpankowski, A. Grama, P. Drineas
%
% Last Update : 10/29/2017

a = zeros(m+1,1);

a(1) = (1/2)*u*(log(u/4)+ 1);

if m>0
    a(2) = (u/4)*(log((u/4)^2)+3);
    for i = 3:m+1
        a(i) = (((-1)^(i-1))*u) / ((i-1)*((i-1)^2-1));
    end
end