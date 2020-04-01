clear

t  = 0:500;
TR = 2;
LowFreq  = sin(2*pi*(1/250)*t);
HighFreq = 3*sin(2*pi*(1/10)*t);
DC = 3; 

Y = DC + HighFreq + LowFreq;
T = numel(Y);

% similar to AFNI
pnum = 1+floor(TR.*T/150);
Yt2 = multpolyfit(repmat(1:T,1,1)',Y',T,pnum);

% matlab easy function -- sanity check
Yt3 = detrend(Y,pnum); % this is useless because only assumes linear in my current matlab version (R2016b)

%Just to check whether xt3 matches the matlab functions
Yt4 = Y - polyval(polyfit(t,Y,pnum),t);

% FSL what does
filter=hp_fsl(T,100,TR);
Yft0 = filter*Y';

SPMdefaultHPCutOff = 128;
% built-in SPM
%K = struct('RT',TR,'row',1:T,'HParam',SPMdefaultHPCutOff);
%Ydct = spm_filter(K,Y'); %DCT is happenening withing the function 
% Similar:
DCT_pnum  = fix(2*(T*TR)/SPMdefaultHPCutOff + 1);
DCT_hpf   = spm_dctmtx(T,DCT_pnum);
dctY = Y' - DCT_hpf*(DCT_hpf'*Y');

figure; hold on; 
plot(t,Y,'LineWidth',1.5,'LineStyle','-.')
plot(t,Yt2,'LineWidth',1.5,'LineStyle','-.','Color','k')
plot(t,Yft0,'LineWidth',1.5,'LineStyle','-.','Color','b')
plot(t,Yt3)
plot(t,Yt4)
plot(t,dctY,'LineWidth',1.5,'LineStyle','-.','Color','r')
legend('raw','poly','hpf','matlab','polyfitmatlabfunc','DCT')
%plot(t,x-y)

mean((HighFreq'-dctY).^2)
mean((HighFreq-Yft0').^2)
mean((HighFreq-Yt2).^2)
