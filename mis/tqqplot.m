function [x,y] = tqqplot(tstat,df)

tstat   = sort(tstat);
np      = 1000;
p       = (0+1/(np-1) : 1/(np-1) : 1-1/(np-1))';
%p(1:5)
x       = tinv(p,df);
y       = quantile(tstat,p);
