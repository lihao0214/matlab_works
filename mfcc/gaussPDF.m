function prob = gaussPDF(Data, Mu, Sigma)
%
% This function computes the Probability Density Function (PDF) of a
% multivariate Gaussian represented by means and covariance matrix.
%
% Author:	Sylvain Calinon, 2009
%			http://programming-by-demonstration.org
%
% Inputs -----------------------------------------------------------------
%   o Data:  D x N array representing N datapoints of D dimensions.
%   o Mu:    D x K array representing the centers of the K GMM components.
%   o Sigma: D x D x K array representing the covariance matrices of the
%            K GMM components.
% Outputs ----------------------------------------------------------------
%   o prob:  1 x N array representing the probabilities for the
%            N datapoints.

[nbVar,nbData] = size(Data);

Data = Data' - repmat(Mu',nbData,1);%X-[u1';u1';...;u1']=[(x1'-u1');(x2'-u1');...;(xN'-u1')]
prob = sum((Data*inv(Sigma)).*Data, 2);%[p1';p2';...;pN']*A=diag([p1'*A;...;pN'*A]*[p1,...,pN])
%prob=diag(Data*inv(Sigma)*Data');%X-[u1';u1';...;u1']=[(x1'-u1');(x2'-u1');...;(xN'-u1')]
prob = exp(-0.5*prob) / sqrt((2*pi)^nbVar * abs(det(Sigma)));%exp(-0.5*(*))
