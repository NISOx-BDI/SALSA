function [x,y] = llpplot(tstat,df,np)

if nargin < 3; np = 5000; end

pval  = 1-tcdf(tstat,df);
spv   = sort(pval);
p     = linspace(1/(np-1),1-(1/np),np); %(1/(np-1) : 1/(np-1) : 1-1/(np-1));
y     = quantile(spv,p);
x     = (1:length(y))/length(y);
x     = quantile(x,p);