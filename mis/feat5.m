ts_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/filtered_func_data.nii.gz';
tcon_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design.con';
dmat_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design_mat.txt';
path2mask='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mask.nii.gz';
parmat='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mc/prefiltered_func_data_mcf.par';
%feat5/featlib.cc 

addpath('/Users/sorooshafyouni/Home/GitClone/FILM2/mis')
addpath('/Users/sorooshafyouni/Home/matlab/spm12')

[Y,ImgStat] = CleanNIFTI_spm(ts_fname,'demean');
Y = Y';
Y = Y - mean(Y);
T=900;

disp('MC params.')
MCp       = load(parmat); 
MCp       = GenMotionParam(MCp,24); 
X         = [load(dmat_fname) MCp];

disp('hpf')
K         = hp_fsl(size(Y,1),100,0.645);    
X         = K*X;    % high pass filter the design
Y         = K*Y;  % high pass filter the data

X = [ones(T,1) X];
tcon      = zeros(1,size(X,2));
tcon(2)   = 1;

tukey_m   = 30; 
tukey_f   = 0; 

ImgStat   = []; 
path2mask = []; 

[cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = feat(Y,X,tcon,tukey_m,tukey_f,ImgStat,path2mask,1,K);

[PSDx,PSDy]   = DrawMeSpectrum(RES,1);
[WPSDx,WPSDy] = DrawMeSpectrum(WRES,1);

%figure; 
hold on; grid on; 
plot(PSDx,mean(PSDy,2))
plot(WPSDx,mean(WPSDy,2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ts_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/filtered_func_data.nii.gz';
% tcon_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design.con';
% dmat_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design_mat.txt';
% path2mask='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mask.nii.gz';
% parmat='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mc/prefiltered_func_data_mcf.par';
% WMSeg='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/reg/func_wmseg.nii.gz';
% %feat5/featlib.cc 
% 
% addpath('/Users/sorooshafyouni/Home/GitClone/FILM2/mis')
% addpath('/Users/sorooshafyouni/Home/matlab/spm12')
% 
% [Y,ImgStat] = CleanNIFTI_spm(ts_fname,'demean');
% Y = Y';
% Y = Y - mean(Y);
% 
% [~,Idx_wm]       = MaskImg(Y',WMSeg,ImgStat); 
% Idx_wm           = Idx_wm{1};
% Ywm              = Y(:,Idx_wm); 
% Idx_gm           = ~Idx_wm;
% Y                = Y(:,Idx_gm); % time series in GM & CSF
% Y                = Y - mean(Y);
% Y                = Y(:,1:2:end); 
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
% X     = K*X;    % high pass filter the design
% Y     = K*Y;  % high pass filter the data
% 
% X = [ones(T,1) X];
% tcon     = zeros(1,size(X,2));
% tcon(2)  = 1;
% 
% tukey_m   = 30; 
% tukey_f   = 0; 
% 
% [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = feat(Y,X,tcon,tukey_m,tukey_f,ImgStat,path2mask,1,[]);
% 
% [PSDx,PSDy]   = DrawMeSpectrum(RES,1);
% [WPSDx,WPSDy] = DrawMeSpectrum(WRES,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = feat(Y,X,tcon,tukey_m,tukey_f,ImgStat,path2mask,badjflag,K)
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
    
    disp('::feat5::')
    if ~exist('badjflag','var');   badjflag = 0;    end;
    
    if ~exist('K','var') || isempty(K)          
        K  = eye(ntp); 
    else
        disp(['feat5:: filter is incorporated into design.'])
    end    
    
    [cbhat,~,~,stat] = myOLS(Y,X,tcon); % not really needed unless for comparison


    se                   = stat.se;
    tv                   = stat.tval;
    zv                   = stat.zval;

    
    pinvX               = pinv(X); 
    R                   = eye(ntp)-X*pinvX; % residual forming matrix 
    RES                 = R*Y;    
    
    if ~badjflag; R     = [];  end
    
    % calc acf, tukey taper, spatially smooth------------------------------
    if isempty(tukey_m); tukey_m = round(sqrt(ntp)); end
    
    R           = R*K; % inject the filter into R
    acf_tukey   = acf_prep(RES,tukey_m,tukey_f,R,ImgStat,path2mask);

    % make the pwfilter----------------------------------------------------
    W_fft       = establish_pwfilter(acf_tukey,ntp);

    % Prewhiten the time series--------------------------------------------
    WY          = prewhiten_timeseries(Y,W_fft);

    % Prewhiten the design
    Xnoint      = X(:,2:end); % put aside the intercept when prewhitening the design
    Xnoint_fft  = fft_model(Xnoint);

    % refit the model to pre-whitened data -------------------------------
    Wcbhat   = zeros(1,nvox);
    WYhat    = zeros(ntp,nvox); 
    WRES     = zeros(ntp,nvox);
    wse      = zeros(nvox,1); 
    wtv      = zeros(nvox,1);
    wzv      = zeros(nvox,1);
    
    disp('feat5:: Refit the prewhitened model.')
    for iv = 1:nvox
        if ~mod(iv,5000); disp(['feat5:: on voxel: ' num2str(iv)]); end; 
        WXnoint = prewhiten_model(Xnoint_fft,W_fft(:,iv),ntp);
        WX      = [ones(ntp,1),WXnoint]; % add back the intercept 
        
        %GLS
        [Wcbhat(iv),WYhat(:,iv),WRES(:,iv),wstat] = myOLS(WY(:,iv),WX,tcon);

        wse(iv)  = wstat.se;
        wtv(iv)  = wstat.tval;
        wzv(iv)  = wstat.zval;    
        
        % BLUSres -- this will take ages. 
        %WBLUSRES(:,iv) = BLUSres(WY(:,iv),WX,(P+1):T);
    end
    %----------------------------------------------------------------------   
    % Re-iterate the whitening with the origina X -- should not be used ---
%     if feat5repeat
%         disp('feat5:: We iterate again...')
%         % We send the prewhitened Y, with the original X in again
%        [~,~,~,~,~,~,Wcbhat,WYhat,WRES,wse,wtv,wzv] = feat5(WY,X,tcon,-2,ImgStat,path2mask,badjflag); 
%     end
    %----------------------------------------------------------------------
    
    disp('feat5:: done.')

end

function acf_tukey = acf_prep(RES,tukey_m,tukey_f,R,ImgStat,path2mask)
% RES should be TxV. T is time, V voxel. 
    
    acfFWHMl = 5;

    [ntp,nvox]  = size(RES);
    [~,~,acv]   = AC_fft(RES,ntp);
    acv         = acv';
    
    % Find an optimal lag if asked-----------------------------------------
    if tukey_m < 0
        acftmp     = acv./acv(1,:);         
        
        if ~isempty(path2mask) || ~isempty(ImgStat)
            
            disp(['feat5:: Autocovariance is smoothed on ' num2str(acfFWHMl) 'mm and taper on ' num2str(tukey_m) ' lag.'])
            acftmp   = ApplyFSLSmoothing(acftmp',acfFWHMl,ImgStat,path2mask)';
        else
            disp('feat5:: No smoothing is done on the ACF.')
        end           
        
        where2stop = FindBreakPoint(acftmp,ntp);
        
%         disp('# of voxel with flat acf.')
%         sum(where2stop==1)
        
        if tukey_m == -1
            
            lp  = 25; 
            hp  = 75; 
            acl = sum(acftmp.^2);
            hidx = acl>prctile(acl,hp);
            lidx = acl<prctile(acl,lp);
            midx = prctile(acl,lp)<=acl & acl<=prctile(acl,hp);
            
            where2stop(midx) = where2stop(midx).*2;
            where2stop(hidx) = where2stop(hidx).*3;
        end
            
        tukey_m    = fix(prctile(where2stop,99.99)); % the max, but avoid outlier
        disp(['feat5:: mean breakpint: ' num2str(mean(where2stop)) ', max:' num2str(max(where2stop)) ', 99th: ' num2str(tukey_m)])
    end
    
    % adjust for bias------------------------------------------------------
    disp('feat5:: adjusting autocovariances: start.')
    if  ~isempty(R)
        if ~ismatrix(R); error('feat5:: adjusting autocovariances: filter should be a matrix'); end; 
        invM = ACFBiasAdjMat(R,ntp,tukey_m);
    else
        invM = eye(tukey_m+1); 
    end

    disp('feat5:: adjusting autocovariances: done.')
    
    acv         = invM*acv(1:tukey_m+1,:); %apply adjustment
    acf         = acv./acv(1,:); % get ACF
    
    % Spatially Smooth ACF-------------------------------------------------
    if ~isempty(path2mask) || ~isempty(ImgStat)
        disp(['feat5:: Autocovariance is smoothed on ' num2str(acfFWHMl) 'mm and taper on ' num2str(tukey_m) ' lag.'])
        acf   = ApplyFSLSmoothing(acf',acfFWHMl,ImgStat,path2mask)';
    else
        disp('feat5:: No smoothing is done on the ACF.')
    end            
    %----------------------------------------------------------------------    
    
    % Tukey Tapering-------------------------------------------------------
    if tukey_m == -1 % this is trouble.
        for i = 1:size(acf,2)
            acfmask                  = zeros(1:ntp); 
            acfmask(1:where2stop(i)) = 1;
            acf(:,i)                 = acf(:,i).*acfmask;
        end
    else
        % Tukey taper
        if tukey_f
            disp(['feat5:: Tukey regularisation.'])
            acf         = acf(1:tukey_m,:);
            lag         = 0:tukey_m-1;
            window      = .5 * (1 + cos(pi .* lag ./ tukey_m));
            acf_tukey   = acf .* repmat(window',1,nvox);
        else
            disp(['feat5:: No Tukey regularisation is set.'])
            acf_tukey   = acf(1:tukey_m,:);
        end
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
    % Does not change the specturm obv.
    %W_fft0 = sqrt(sum(W_fft(2:end,:).^2))./z_pad;
    %W_fft = W_fft./W_fft0;

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

function invM_biasred = ACFBiasAdjMat(R,ntp,ARO)
% Bias adjustment for ACF of residuals. 
% 
% This is a bit tricky here. Note that ResidFormingMat is symmetric. 
% So, R'*Di*R*Dj == R'*Di*R*Dj
% And therefore, we can have R = R*K; if we wanted to inject the K filter into the R. 
% 
% SA, Ox, 2020

% fmristat implementation ------------------------------------------------
% of implementation. 
%     M_biasred   = zeros(ARO+1);
%     for i=1:(ARO+1)
%         Di                  = (diag(ones(1,ntp-i+1),i-1)+diag(ones(1,ntp-i+1),-i+1))/(1+(i==1));
%         for j=1:(ARO+1)
%            Dj               = (diag(ones(1,ntp-j+1),j-1)+diag(ones(1,ntp-j+1),-j+1))/(1+(j==1));
%            M_biasred(i,j)   = trace(ResidFormingMat'*Di*ResidFormingMat*Dj)/(1+(i>1));
%         end
%     end
%     invM_biasred = inv(M_biasred);
% ------------------------------------------------------------------------

% Appendix A & then mofidied for K from the MS Notes ---------------------
    M_biasred   = zeros(ARO+1);
    for l = 1:(ARO+1)
        Dl = diag(ones(1,ntp-l+1),l-1); % upper triangle
        for j = 1:(ARO+1)
           DjDjt            = (diag(ones(1,ntp-j+1),j-1)+diag(ones(1,ntp-j+1),-j+1))/(1+(j==1));
           M_biasred(l,j)   = trace(R'*Dl*R*DjDjt);
        end
    end
    invM_biasred = inv(M_biasred);
% ------------------------------------------------------------------------

end

function where2stop = FindBreakPoint(acf,T)
% this finds the breaking points for shrinking the AC. 
% Nothing serious, just might help with speed...
% SA, Ox, 2018-2020
    if ~sum(ismember(size(acf),T)); error('There is something wrong!'); end
    if size(acf,2) ~= T; acf = acf'; end
    
    bnd        = (sqrt(2)*erfinv(0.95))./sqrt(T); %CI of ACF
    isit       = abs(acf)>bnd;

    where2stop = zeros(1,size(acf,1));
    for i = 1:size(acf,1)
        where2stop(i) = find(~isit(i,:),1)-1;
    end
end

