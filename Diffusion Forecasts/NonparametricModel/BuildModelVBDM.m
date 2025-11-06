function [basis,ForwardOp,peqoversample,peq,qest,normalizer,X,epsilon] = BuildModelVBDM(x,k,k2,nvars,delays,shifts,dim)
%%% Written by Tyrus Berry of George Mason University

%%% Inputs: x       - T x N time series of length T with dimension N
%%%         k       - number of nearest neighbors to use
%%%         k2      - nearest neighbors to use in pre-bandwidth
%%%         nvars   - number of eigenfunctions to find
%%%         delays  - number of delays to use
%%%         shifts  - number of shifts to use in the shift operator
%%%         dim     - (optional) intrinsic dimension of the attractor

%%% Outputs: basis     - eigenfunctions of gradient flow on attractor
%%%          ForwardOp - nvars x nvars forecast model in eigencoords.
%%%          peq       - invariant measure
%%%          qest      - sampling measure
%%%          normalizer- tool to normalize coefficients of a density
%%%          X         - delay embedded data
%%%          epsilon   - diffusion maps parameter (auto fit)

T=size(x,1);
N=size(x,2);
L = length(delays);
maxd = max(delays);

X=zeros(T-maxd+1,N*L);

for i=1:L
    X(:,(i-1)*N + (1:N)) = x(delays(i):T-maxd+delays(i),:);
end

operator = 2; %%% build the generator for the gradient flow system

if (nargin < 7)
    [basis,eigs,epsilon,peqoversample,peq,qest] = VBDM(X,k,k2,nvars,operator);
else 
    [basis,eigs,epsilon,peqoversample,peq,qest] = VBDM(X,k,k2,nvars,operator,dim);
end

lambdas = -log(diag(eigs));

normalizer = mean(basis.*repmat(peq./qest,1,nvars));

T=T-maxd;

ForwardOp = (basis(1+shifts:T,:))'*(basis(1:T-shifts,:).*repmat(peq(1:T-shifts)./qest(1:T-shifts),1,nvars))/(T-shifts);







