function H = TaylorEntropy(R,m,p,u)
% This function computes the Entropy of the density matrix R, using the
% Taylor approximation and the **Gaussian** trace estimator. Method
% appeared in [1].
%
% Inputs:
%
% R       : n x n input matrix (density matrix)
% m       : number of Taylor Terms to keep
% p       : number of Gaussian vectors
% u       : estimation to the largest singular value / eigenvalue
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
% -- Last Update 10/6/2020

n     = size(R,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u = min([6*(powermethod(R,ceil(log(n)),ceil(log(1/delta)))),1]); 
% estimate the largest eigenvalue
% u     = min([6*eigs(R,1);1]);
% u     = 1;
% msg   = sprintf('\t Estimation of largest eigenvalue = % f',u); disp(msg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C     = eye(n) - (1/u).*R;
f1    = log(1/u);

G     = randn(n,p); % Gaussian vectors 
GAMMA = zeros(p,1);

for i = 1:p	
    ls             = G(:,i)'*R; % this is the left side of the taylor estimation
    v              = C*G(:,i);  % (I-inv(u)R)*g_i
    temp           = ls*v;      % 1st term of the truncated taylor series
    for k = 2:m       
	v          = C*v;       % nominator of the second (e.t.c.) term of the truncated series 
        temp       = temp + (ls*v)/k; % taylor series
    end
    GAMMA(i)= temp;
end

f2 = mean(GAMMA);
fprintf('\t Trace Estimator = %f\n',f2);

H = f1 + f2;
fprintf('\t Entropy = %f\n', H);
