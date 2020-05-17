% ts_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/filtered_func_data.nii.gz';
% tcon_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design.con';
% dmat_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design_mat.txt';
% path2mask='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mask.nii.gz';
% parmat='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mc/prefiltered_func_data_mcf.par';
% %feat5/featlib.cc 
% 
% [Y,ImgStat] = CleanNIFTI_spm(ts_fname,'demean');
% Y = Y';
% Y = Y - mean(Y);
% 
% disp('MC params.')
% MCp      = load(parmat); 
% MCp      = GenMotionParam(MCp,24); 
% X        = [ones(T,1) load(dmat_fname) MCp];
% 
% tcon     = zeros(1,size(X,2));
% tcon(2)  = 1;
% 
% addpath('/Users/sorooshafyouni/Home/GitClone/FILM2/mis')
% 
% disp('hpf')
% hp_ff = hp_fsl(size(Y,1),100,0.645);    
% X     = hp_ff*X;    % high pass filter the design
% Y     = hp_ff*Y;  % high pass filter the data
% 
% [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = feat(Y,X,tcon,[],ImgStat,path2mask,1);
% 
% [PSDx,PSDy]   = DrawMeSpectrum(RES,1);
% [WPSDx,WPSDy] = DrawMeSpectrum(WRES,1);
% 
% figure; hold on; grid on; 
% plot(PSDx,mean(PSDy,2))
% plot(WPSDx,mean(WPSDy,2))


function [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = feat5(Y,X,tcon,tukey_m,ImgStat,path2mask,badjflag)
% Y      : TxV
% X      : TxEV
% tcon   : 1xEV
% tukey_m: int
%
% Y should gone through detrending & smoothing prior to this function.
%
% ${FSLDIR}/src/feat5/prewhiten.cc
% ${FSLDIR}/src/feat5/featlib.cc
%
% SA, Ox, 2020

    if ~exist('badjflag','var'); badjflag = 0; end;

    ntp                   = size(Y,1); 
    nvox                  = size(Y,2); 
    [cbhat,~,~,stat] = myOLS(Y,X,tcon); % not really needed unless for comparison


    se                   = stat.se;
    tv                   = stat.tval;
    zv                   = stat.zval;

    
    pinvX               = pinv(X); 
    ResidFormingMat     = eye(ntp)-X*pinvX; % residual forming matrix 
    RES                 = ResidFormingMat*Y;    
    
    if ~badjflag; ResidFormingMat = [];  end
    
    % calc acf & tukey taper it
    if isempty(tukey_m); tukey_m = round(sqrt(ntp)); end
    
    acf_tukey   = acf_prep(RES,tukey_m,ResidFormingMat);

    % smooth ACF
    FWHMl       = 5; 
    acf_tukey   = ApplyFSLSmoothing(acf_tukey',FWHMl,ImgStat,path2mask)';
    % make the pwfilter
    W_fft       = establish_pwfilter(acf_tukey,ntp);

    % Prewhiten the time series
    WY          = prewhiten_timeseries(Y,W_fft);

    % Prewhiten the design
    X_fft       = fft_model(X);

    % refit the model to pre-whitened data
    Wcbhat = zeros(1,nvox);
    WYhat  = zeros(ntp,nvox); 
    WRES   = zeros(ntp,nvox);
    wse    = zeros(nvox,1); 
    wtv    = zeros(nvox,1);
    wzv    = zeros(nvox,1);
    for iv = 1:nvox
        if ~mod(iv,5000); disp(['on voxel: ' num2str(iv)]); end; 
        WX                        = prewhiten_model(X_fft,W_fft(:,iv),ntp);
        [Wcbhat(iv),WYhat(:,iv),WRES(:,iv),wstat] = myOLS(WY(:,iv),WX,tcon);

        wse(iv)  = wstat.se;
        wtv(iv)  = wstat.tval;
        wzv(iv)  = wstat.zval;    

    end

end

function acf_tukey = acf_prep(RES,tukey_m,ResidFormingMat)
% RES should be TxV
    
    [ntp,nvox]  = size(RES);
    [~,~,acv]   = AC_fft(RES,ntp);
    acv         = acv';
    
    % adjust for bias
    if nargin==3 || isempty(ResidFormingMat)
        BiasAdj = ACFBiasAdjMat(ResidFormingMat,ntp,tukey_m);
        disp('feat:: adjusting autocovariances.')
    elseif nargin<3
        BiasAdj = eye(tukey_m+1); 
    end
    
    
    acv         = BiasAdj*acv(1:tukey_m+1,:); %apply adjustment
    acf         = acv./acv(1,:); % get ACF
    
    
    % Tukey taper
    acf         = acf(1:tukey_m,:);
    lag         = 0:tukey_m-1;
    window      = .5 * (1 + cos(pi .* lag ./ tukey_m));
    acf_tukey   = acf .* repmat(window',1,nvox);
    
end

function W_fft = establish_pwfilter(acf,ntp)
% acf : MxV matrix. M is the tukey M
    nvox                              = size(acf,2); 
    z_pad                             = 2.^ceil(log(ntp)/log(2));
    tukey_m                           = size(acf,1);
    acf_kernel                        = zeros(z_pad, nvox);
    acf_kernel(1:tukey_m,:)           = acf;
    acf_kernel((end-tukey_m)+2:end,:) = flip(acf(2:end,:));
    
    acf_fft                           = real(fft(acf_kernel,[],1));
    W_fft                             = zeros(z_pad, nvox);
    W_fft(2:end,:)                    = 1 ./ sqrt(abs(acf_fft(2:end,:)));
    
    % feat5/featlib.cc  line 127:
    % But why? If we do this, we'll lose the variance in the original signal
    %W_fft0 = sqrt(sum(W_fft(2:end,:).^2))./w_pad;
    %W_fft = W_fft./W_fft0;

end

function WY = prewhiten_timeseries(Y,W_fft)
% Y is TxV
    ntp   = size(Y,1);
    z_pad = 2.^ceil(log(ntp)/log(2));
    Y_fft = fft(Y,z_pad,1);
    WY    = real(ifft((W_fft.*real(Y_fft)+W_fft.*imag(Y_fft)*1j),[],1));
    WY    = WY(1:ntp,:);
end


function X_fft = fft_model(X)
% X is TxEV
    [ntp,nev]  = size(X);
    z_pad      = 2.^ceil(log(ntp)/log(2));
    X_fft      = zeros(z_pad,nev);
    for ev_cnt = 1:nev
        X_fft(:,ev_cnt) = fft(X(:,ev_cnt),z_pad,1);
    end
end

function WX = prewhiten_model(X_fft,W_fft_v,ntp)
    [z_pad,nev]  = size(X_fft);
    WX           = zeros(ntp,nev);
    for ev_cnt = 1:nev
        X_fft_ev = X_fft(:,ev_cnt);
        WX0          = real(ifft((W_fft_v.*real(X_fft_ev)+W_fft_v.*imag(X_fft_ev)*1j),z_pad,1));
        WX(:,ev_cnt) = WX0(1:ntp,:);
    end
end

function invM_biasred = ACFBiasAdjMat(ResidFormingMat,ntp,ARO)
% Bias adjustment for ACF of residuals. 

    M_biasred   = zeros(ARO+1);
    for i=1:(ARO+1)
        Di                  = (diag(ones(1,ntp-i+1),i-1)+diag(ones(1,ntp-i+1),-i+1))/(1+(i==1));
        for j=1:(ARO+1)
           Dj               = (diag(ones(1,ntp-j+1),j-1)+diag(ones(1,ntp-j+1),-j+1))/(1+(j==1));
           M_biasred(i,j)   = trace(ResidFormingMat*Di*ResidFormingMat*Dj)/(1+(i>1));
        end
    end
    invM_biasred = inv(M_biasred);
    
end


