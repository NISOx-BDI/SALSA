clear

% T    = 500;
% TR   = 3; 
% Fs   = 1/TR;
% nRlz = 5000; 
% %simulate data using filter
% a    = [1 -0.5];%ARcoeffs
% b    = 1 ;%MAcoeffs
% 
% for i = 1:nRlz
%     e = randn(T,1); % generate gaussian white noise
%     Yt(:,i) = filter(b,a,e); % generate y
% end
% 
% %Yt = Yt./std(Yt);
% Yt = Yt - mean(Yt);

Path2Mat = ['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/ROCKLAND/sub-A00028858/ses-DS2/sub-A00028858_ses-DS2_task-rest_acq-645_bold_mpp/dY.mat'];
img = load(Path2Mat);
Yt = img.dY; 
TR = img.TR;

Yt = Yt(:,1:20:end); 

[T,V] = size(Yt) 
Yt = Yt-repmat(mean(Yt),900,1);

disp('Prewhiten with SPD')
ts0 = PrewhitenSPD(Yt,20,1/TR);

[xx,yy] = DrawMeSpectrum(ts0,TR); 


disp('Prewhiten with ACF')
[~,~,Ytacov]  = AC_fft(Yt,T); % Autocovariance; VxT
Ytacov        = Ytacov'; %TxV

for vi = 1:V
    if ~mod(vi,100); disp(['on voxel ' num2str(vi)]); end; 
    [sqrtmVhalf,spdflag] = ACF_ResPWm(Ytacov(:,vi),30,[],1);
    pwYt(:,vi) = sqrtmVhalf*Yt(:,vi);
end


[pwxx,pwyy] = DrawMeSpectrum(pwYt,TR); 

figure; 
hold on; grid on; 
plot(xx,mean(yy,2),'linewidth',1.3)
plot(pwxx,mean(pwyy,2),'LineWidth',1.3)
plot([0 0.7],[1 1],'color',[.5 .5 .5],'LineWidth',1.3)
legend({'Spectral(AR=20)','ACF(Lag=30)'})
ylabel('Power')
xlabel('Freq')