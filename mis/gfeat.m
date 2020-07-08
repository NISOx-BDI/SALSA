% ts_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/filtered_func_data.nii.gz';
% tcon_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design.con';
% dmat_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design_mat.txt';
% path2mask='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mask.nii.gz';
% parmat='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mc/prefiltered_func_data_mcf.par';
% WMSeg='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/reg/func_wmseg.nii.gz';
% %feat5/featlib.cc 
% 
% addpath('/Users/sorooshafyouni/Home/GitClone/FILM2/mis')
% addpath('/Users/sorooshafyouni/Home/GitClone/FILM2/utils/Trend')
% addpath('/Users/sorooshafyouni/Home/matlab/spm12')
% 
% [Y,ImgStat] = CleanNIFTI_spm(ts_fname,'demean');
% Y = Y';
% Y = Y - mean(Y);
% 
% % [~,Idx_wm]       = MaskImg(Y',WMSeg,ImgStat); 
% % Idx_wm           = Idx_wm{1};
% % Ywm              = Y(:,Idx_wm); 
% % Idx_gm           = ~Idx_wm;
% % Y                = Y(:,Idx_gm); % time series in GM & CSF
% % Y                = Y - mean(Y);
% 
% %%Y                = Y(:,1:20:end); 
% %%size(Y)
% 
% T=900;
% TR=0.645;
% disp('MC params.')
% MCp      = load(parmat); 
% MCp      = GenMotionParam(MCp,24); 
% X        = [load(dmat_fname) MCp];
% 
% disp('hpf')
% K = hp_fsl(size(Y,1),100,0.645);    
% X = K*X;    % high pass filter the design
% Y = K*Y;  % high pass filter the data
% 
% X        = [ones(T,1) X];
% tcon     = zeros(1,size(X,2));
% tcon(2)  = 1;
% 
% aclageval = 0; 
% tukey_m   = sqrt(T); 
% tukey_f   = 1;
% poolflag  = 0; 
% [WY,WX,cbhat,RES,ostat,se,tv,zv,Wcbhat,WRES,wse,wtv,wzv] = gfeat5(Y,X,TR,tcon,tukey_m,tukey_f,aclageval,1,K,poolflag);
% 
% [PSDx,PSDy]   = DrawMeSpectrum(RES,1);
% [WPSDx,WPSDy] = DrawMeSpectrum(WRES,1);
% 
% %figure; 
% hold on; grid on; 
% %plot(PSDx,mean(PSDy,2))
% plot(WPSDx,mean(WPSDy,2))

function [WY,WX,cbhat,RES,ostat,se,tv,zv,Wcbhat,WRES,wse,wtv,wzv] = gfeat(Y,X,TR,tcon,tukey_m,tukey_f,aclageval,badjflag,K,poolflag)
% Y      : TxV
% X      : TxEV. Always always intercept is the first column
% tcon   : 1xEV
% tukey_m: int
%
% Y should gone through detrending & smoothing prior to this function.
%
% ${FSLDIR}/src/feat5/prewhiten.cc
% ${FSLDIR}/src/feat5/featlib.cc
%
% SA, Ox, 2020

    [ntp,nvox]  = size(Y); 
    
    disp('::gfeat::')
    if ~exist('badjflag','var');   badjflag = 0;    end;
       
    if ~exist('aclageval','var') || isempty(aclageval) && tukey_m ~=-2
        aclageval = tukey_m.^2; 
    elseif ~exist('aclageval','var') || isempty(aclageval) && tukey_m ==-2
        error('gfeat:: Specify aclageval.')
    end; 

    if ~exist('K','var') || isempty(K)          
        K  = eye(ntp); 
    else
        disp(['gfeat:: filter is incorporated into design.'])
    end     

    disp(['gfeat:: #lags used to estimate acf: ' num2str(tukey_m)])
    disp(['gfeat:: #lags used to adjust acf:   ' num2str(aclageval)])    

    
    pinvX              = pinv(X); 
    R                  = eye(ntp)-X*pinvX; % residual forming matrix 
    RES                = R*Y;    
    R                  = R*K; % inject the filter into R
   %%%%% REMOVE THIS AFTER EVALUATION
    [cbhat,~,~,stat] = myOLS(Y,X,tcon);
    se                 = stat.se;
    tv                 = stat.tval;
    zv                 = stat.zval;
    ostat.df           = stat.df; 
    clear stat
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % calc acf & tukey taper it
    if isempty(tukey_m); tukey_m = round(sqrt(ntp)); end
    if ~badjflag;        R = [];                     end
    if ~exist('poolflag','var') || isempty(poolflag); poolflag = 1; end 
    
    acf_tukey   = g_acf_prep(RES,TR,tukey_m,tukey_f,R,aclageval,poolflag);

    % make the pwfilter
    W_fft       = establish_pwfilter(acf_tukey,ntp);
    W_fft       = median(W_fft,2); % median across voxels, just to avoid outliers. 
    
    % Prewhiten the time series
    WY          = prewhiten_timeseries(Y,W_fft);

    % Prewhiten the design
    Xnoint      = X(:,2:end); % put aside the intercept when prewhitening the design
    Xnoint_fft  = fft_model(Xnoint);

    WXnoint     = prewhiten_model(Xnoint_fft,W_fft,ntp);
    WX          = [ones(ntp,1),WXnoint]; % add back the intercept 

    %GLS
    [Wcbhat,~,WRES,wstat] = myOLS(WY,WX,tcon);
    wse                   = wstat.se;
    wtv                   = wstat.tval;
    wzv                   = wstat.zval;    
    clear wstat
    
    % refit the model to pre-whitened data -------------------------------
%     Wcbhat   = zeros(1,nvox);
%     WYhat    = zeros(ntp,nvox); 
%     WRES     = zeros(ntp,nvox);
%     wse      = zeros(nvox,1); 
%     wtv      = zeros(nvox,1);
%     wzv      = zeros(nvox,1);
%     
%     disp('fastfeat:: Refit the prewhitened model.')
%     for iv = 1:nvox
%         if ~mod(iv,5000); disp(['fastfeat:: on voxel: ' num2str(iv)]); end; 
%         WXnoint = prewhiten_model(Xnoint_fft,W_fft(:,iv),ntp);
%         WX      = [ones(ntp,1),WXnoint]; % add back the intercept 
%         
%         %GLS
%         [Wcbhat(iv),WYhat(:,iv),WRES(:,iv),wstat] = myOLS(WY(:,iv),WX,tcon);
% 
%         wse(iv)  = wstat.se;
%         wtv(iv)  = wstat.tval;
%         wzv(iv)  = wstat.zval;    
%         
%         % BLUSres -- this will take ages. 
%         %WBLUSRES(:,iv) = BLUSres(WY(:,iv),WX,(P+1):T);
%     end
       
    disp('gfeat:: done.')

end

function acf_tukey = g_acf_prep(RES,TR,tukey_m,tukey_f,R,aclageval,poolflag)
% RES should be TxV. T is time, V voxel. 
    

    [ntp,nvox]      = size(RES);
    [~,acfCI,acv]   = AC_fft(RES,ntp);
    acv             = acv';
    
    if poolflag
        disp(['gfeat:: pooling by autocorrelation.'])
        %acl          = sum(acf(1,1:fix(ntp/4)).^2); % Anderson's suugestion of ignoring beyond ntps/4
        acf1        = acv(1+1,:)./acv(1,:); %acf(1) 
        jidx        = find(abs(acf1)>acfCI(2));  % only if a voxel exceeds the CI    
        nvox        = numel(jidx);
        acv         = acv(:,jidx);
    else
       disp('gfeat:: No pooling is done.') 
    end
    
    if ~nvox
        error('gfeat:: Something is wrong, there should be at least some voxels significant to the overall design'); 
    else
        disp(['gfeat:: Number of pooled voxels: ' num2str(nvox)])
    end    
    
    if tukey_m < 0
        %flag2      = 1; 
        acftmp     = acv./acv(1,:);        
     
        where2stop = FindBreakPoint(acftmp,ntp);
        clear acftmp
   
        tukey_m    = fix(prctile(where2stop,99.999)); % the max, but avoid outlier
        disp(['gfeat:: mean breakpint: ' num2str(mean(where2stop)) ', max:' num2str(max(where2stop)) ', 99th: ' num2str(tukey_m)])
        where2stop(where2stop>tukey_m) = tukey_m; % set everything above tukey_m to tukey_m
    end
    
    maxlag      = max(aclageval,tukey_m);
    acv         = acv(1:maxlag+1,:);    
        
    % adjust for bias
    disp('gfeat:: adjusting autocovariances: start.')
    if  ~isempty(R)
        if ~ismatrix(R); error('gfeat:: adjusting autocovariances: filter should be a matrix'); end; 
        if aclageval
            M = ACFBiasAdjMat(R,ntp,aclageval);
        else
            M = ACFBiasAdjMat(R,ntp,tukey_m);
        end
    else
        M = eye(aclageval+1); 
    end
    
    
    if aclageval 
        Bfull = ACFbasis(ntp-1,TR); 
        Best  = Bfull(1:tukey_m+1,:);
        Beval = Bfull(1:aclageval+1,:);

        g     = (M*Beval)\acv(1:aclageval+1,:);
        acf1   = Best*g;
        acf1   = acf1./acf1(1,:); 
    else
        acv  = pinv(M)*acv(1:tukey_m+1,:); %apply adjustment
        acf1  = acv./acv(1,:); % get ACF        
    end
    
    disp('gfeat:: adjusting autocovariances: done.')

    % Tukey taper
    if tukey_f
        disp(['gfeat:: Tukey regularisation.'])
        acf1                        = acf1(1:tukey_m,:);
        lag                         = 0:tukey_m-1;
        window                      = .5 * (1 + cos(pi .* lag ./ tukey_m));
        acf_tukey                   = acf1 .* repmat(window',1,nvox);
    else
        disp(['gfeat:: No Tukey regularisation is set.'])
        acf_tukey                   = acf1(1:tukey_m,:);
    end

      
end

function W_fft = establish_pwfilter(acf,ntp)
% acf : MxV matrix. M is the tukey M. V is voxel
    nvox                              = size(acf,2); 
    z_pad                             = 2.^ceil(log(ntp)/log(2)); % that is 2.^nextpow2 
    tukey_m                           = size(acf,1);
    acf_kernel                        = zeros(z_pad, nvox);
    acf_kernel(1:tukey_m,:)           = acf;
    acf_kernel((end-tukey_m)+2:end,:) = flip(acf(2:end,:));
    
    acf_fft                           = real(fft(acf_kernel,[],1));
    W_fft                             = zeros(z_pad, nvox);
    W_fft(2:end,:)                    = 1 ./ sqrt(abs(acf_fft(2:end,:)));
    
    % feat5/featlib.cc  line 127:
    % But why? If we do this, we'll lose the variance in the original signal
    W_fft0 = sqrt(sum(W_fft(2:end,:).^2))./z_pad;
    W_fft = W_fft./W_fft0;

end

function WY = prewhiten_timeseries(Y,W_fft)
% Y is TxV. T is time, V is voxel.
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

 function M = ACFBiasAdjMat(R,ntp,ARO)
% Bias adjustment for ACF of residuals. 
% R   : residual forming matrix
% ntp : # of data points
% ARO : AR order 
%
% --
% This is a bit tricky here. Note that ResidFormingMat is symmetric. 
% So, R'*Di*R*Dj == R'*Di*R*Dj
% And therefore, we can have R = R*K; if we wanted to inject the K filter into the R. 
% 
% SA, Ox, 2020

% Appendix A & then mofidied for K from the MS Notes ---------------------
    disp(['gfeat:: BiasAdjMat::'])
    M   = zeros(ARO+1);
    for l = 1:(ARO+1)
        if ~mod(l,50); disp(['gfeat:: BiasAdjMat:: on lag: ' num2str(l)]); end; 
        Dl = diag(ones(1,ntp-l+1),l-1); % upper triangle
        for j = 1:(ARO+1)
           if j < l; continue; end; % M is symm so, save time!
           DjDjt    = (diag(ones(1,ntp-j+1),j-1)+diag(ones(1,ntp-j+1),-j+1))/(1+(j==1));
           M(l,j)   = trace(R'*Dl*R*DjDjt);
        end
    end
    
    dM = M.*~eye(size(M)); % set the diag to zero
    M  = M + dM'; % add back the lower triangle
    
% ------------------------------------------------------------------------

 end

function where2stop = FindBreakPoint(acf,T)
% this finds the breaking points for shrinking the AC. 
% Nothing serious, just might help with speed...
% SA, Ox, 2018-2020
    if ~sum(ismember(size(acf),T)); error('gfeat:: where2stop:: There is something wrong!'); end
    if size(acf,2) ~= T; acf = acf'; end
    
    bnd        = (sqrt(2)*erfinv(0.95))./sqrt(T); %CI of ACF
    isit       = abs(acf)>bnd;

    where2stop = zeros(1,size(acf,1));
    for i = 1:size(acf,1)
        where2stop(i) = find(~isit(i,:),1)-1;
    end
end

 function B = ACFbasis(MaxLag,TR)
 % matrix of Mxlag x basis
 % TEN, 2020
    t = [(0:MaxLag)*TR]';
    d = 2.^(floor(log2(TR/4)):log2(64));
    d(7:end) = []; % as per spm_Ce.m code.
    j = [0 1];%[0 1]; % or	[0 1 2], but I don't think j=2 buys us much
    J = repmat(j,length(t),length(d));
    T = repmat(t,1,length(j)*length(d));
    D = repmat(repelem(d,length(j)),length(t),1);
    B = T.^J.*exp(-T./D);
    B = B./std(B);
 end
