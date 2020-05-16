

clear
nRlz = 10000; 
T    = 1000-1; 

% nts = 114;
% ts  = randn(T,nts)+randn(T,1);
% idx = find(tril(ones(nts),-1));
% 
% ts = ts-mean(ts);
% ts = ts./std(ts);
% 
% [xx,yy] = DrawMeSpectrum(ts,1);
% 
% plot(xx,mean(yy,2))
% 
% rts = corr(ts);
% zts = atanh(rts(idx)).*sqrt(T-3);
% 
% sts = ThetaRnd(ts,T);
% rsts1 = corr(sts);
% zsts1 = atanh(rsts1(idx)).*sqrt(T-3);
% 
% 
% sts = ThetaRnd(ts,T,0);
% rsts0 = corr(sts);
% zsts0 = atanh(rsts0(idx)).*sqrt(T-3);
% 
% % addpath('/Users/sorooshafyouni/Home/GitClone/xDF')
% % [A,B]  = xDF(sts,T,'taper','tukey',sqrt(T));
% % zsts3 = B.z(idx);
% 
% sts00 = phaseRandomize(ts(1:T,:));
% rsts4 = corr(sts00);
% zsts4 = atanh(rsts4(idx)).*sqrt(T-3);


a    = [1 -0];%ARcoeffs
b    = [1 0] ;%MAcoeffs

rho = 0;
for i = 1:nRlz
    % i.i.d cross correlated time series
    %ts = mvnrnd([0 0],[1 rho;rho 1],T);
    
    % autocorrelated decorrelated time series
    R = [1 rho;rho 1];
    R = chol(R);
    
    ts = filter(b,a,randn(T,2));
    
    ts = ts*R';
    
    cc = atanh(corr(ts)).*sqrt(T-3);
    z_naive(i) = cc(1,2);
    
    ts0 = ThetaRnd(ts,T);
    cc0 = atanh(corr(ts0)).*sqrt(T-3);
    z_naive0(i) = cc0(1,2);
    
    ts1 = ThetaRnd(ts,T,0);
    cc1 = atanh(corr(ts1)).*sqrt(T-3);
    z_naive1(i) = cc1(1,2);        
    
end


figure; hold on; 

[x,y] = HistLine_pdf(z_naive,50);
plot(x,y,'linewidth',1.7,'color','b')
[x,y] = HistLine_pdf(z_naive0,50);
plot(x,y,'linewidth',1.7,'color','r')
[x,y] = HistLine_pdf(z_naive1,50);
plot(x,y,'linewidth',1.7,'color','k')

XX = -10:(1+1)/100:10;
plot(XX,normpdf(XX),'color',[0.5 0.5 0.5 0.5],'linewidth',1.7)

figure; hold on; 
[x,y] = nqqplot(z_naive);  plot(x,y)
[x,y] = nqqplot(z_naive0); plot(x,y)
[x,y] = nqqplot(z_naive1); plot(x,y)
abline(0,1)

% [x,y] = HistLine_pdf(zsts0,20);
% plot(x,y)

% [x,y] = HistLine_pdf(zsts3,20);
% plot(x,y)
% [x,y] = HistLine_pdf(zsts4,20);
% plot(x,y)
% [x,y] = HistLine_pdf(zsts22,20);
% plot(x,y,'linewidth',1.7,'color','r')

% [x,y] = HistLine_pdf(zstsRnd,100);
% plot(x,y)

%legend({'Raw corr','Phase Rand. with Corr','Phase Rand. without Corr','xDF','N(0,1)'})

% histogram(rts,100,'Normalization','probability')
% histogram(rsts1,100,'Normalization','probability')
% histogram(rsts0,100,'Normalization','probability')
% legend({'Raw corr','Phase Rand. with Corr','Phase Rand. without Corr'})