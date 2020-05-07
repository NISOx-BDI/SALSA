function [x,y] = nqqplot(tstat)

tstat   = sort(tstat);
np      = 5000;
p       = (0+1/(np-1) : 1/(np-1) : 1-1/(np-1))';
%p(1:5)
x       = norminv(p);
y       = quantile(tstat,p);
