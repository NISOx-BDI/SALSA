clear

t  = 0:300;
TR = 1;
x = 3 + 3*sin(2*pi*(1/10)*t) + sin(2*pi*(1/250)*t);
D = numel(x);

% similar to AFNI
pnum = 1 + round(D/150);
xt2 = detrendnonlin(x,pnum);

% FSL what does
filter=hp_fsl(D,100,TR);
xft2 = filter*x';

figure; hold on; 
plot(t,x)
plot(t,xt2)
plot(t,xft2)
legend('raw','poly','hpf')
%plot(t,x-y)