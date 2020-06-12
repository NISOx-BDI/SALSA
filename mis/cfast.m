clear

ts_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/filtered_func_data.nii.gz';
tcon_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design.con';
dmat_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design_mat.txt';
path2mask='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mask.nii.gz';
parmat='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mc/prefiltered_func_data_mcf.par';
WMSeg='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/reg/func_wmseg.nii.gz';
%feat5/featlib.cc 

addpath('/Users/sorooshafyouni/Home/GitClone/FILM2/mis')
addpath('/Users/sorooshafyouni/Home/GitClone/FILM2/utils/Trend')
addpath('/Users/sorooshafyouni/Home/matlab/spm12')

[Y,ImgStat] = CleanNIFTI_spm(ts_fname,'demean');
Y = Y';
Y = Y - mean(Y);
T=900;
TR = 0.645; 

disp('MC params.')
MCp      = load(parmat); 
MCp      = GenMotionParam(MCp,24); 
X        = [load(dmat_fname) MCp];

X        = X-mean(X); 

disp('hpf')
K = hp_fsl(size(Y,1),100,TR);    
X     = K*X;    % high pass filter the design
Y     = K*Y;  % high pass filter the data

X = [ones(T,1) X];
tcon     = zeros(1,size(X,2));
tcon(2)  = 1;

% [cbhat_gm,RES_gm,stat_gm,se_gm,tv_gm,zv_gm,Wcbhat_gm,WYhat_gm,WRES_gm,wse_gm,wtv_gm,wzv_gm] = nfast5(Ygm,X,0.645,tcon); 
% [PSDx_gm,PSDy_gm]   = DrawMeSpectrum(RES_gm,1);
% [WPSDx_gm,WPSDy_gm] = DrawMeSpectrum(WRES_gm,1);
% [cbhat_wm,RES_wm,stat_wm,se_wm,tv_wm,zv_wm,Wcbhat_wm,WYhat_wm,WRES_wm,wse_wm,wtv_wm,wzv_wm] = nfast5(Ywm,X,0.645,tcon); 
% [PSDx_wm,PSDy_wm]   = DrawMeSpectrum(RES_wm,1);
% [WPSDx_wm,WPSDy_wm] = DrawMeSpectrum(WRES_wm,1);
% figure; hold on; grid on; 
% plot(PSDx_wm,mean(PSDy_wm,2))
% plot(WPSDx_wm,mean(WPSDy_wm,2))

[cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = cfast5(Y,X,TR,tcon,ImgStat,WMSeg);

[PSDx,PSDy]   = DrawMeSpectrum(RES,1);
[WPSDx,WPSDy] = DrawMeSpectrum(WRES,1);

figure; hold on; grid on;
plot(PSDx,mean(PSDy,2))
plot(WPSDx,mean(WPSDy,2))



function [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = cfast5(Y,X,TR,tcon,ImgStat,WMSeg,pmethod)

    if ~exist('pmethod','var');   pmethod = 0;       end; 
    
    [ntp,nvox]          = size(Y);
    
    [cbhat,~,RES,stat]  = myOLS(Y,X,tcon);
    se                  = stat.se;
    tv                  = stat.tval;
    zv                  = stat.zval;
    [~,Idx_wm]          = MaskImg(Y',WMSeg,ImgStat); % Time series in WM 
    Idx_wm              = Idx_wm{1};
    Ywm                 = Y(:,Idx_wm); 
    Idx_gm              = ~Idx_wm;
    Ygm                 = Y(:,Idx_gm); % time series in GM & CSF

    %size(Ygm), size(Ywm)
    
    disp('cfast:: apply gFAST on grey matter')
    chunksize_gm = 10000; 
    [~,~,~,~,~,~, Wcbhat_gm,WYhat_gm,WRES_gm,~,wse_gm,wtv_gm,wzv_gm] = nfast(Ygm,X,TR,tcon,pmethod,chunksize_gm);
    
    chunksize_wm = 10000; 
    disp('cfast:: apply feat5 on the white matter.')
    [~,~,~,~,~,~, Wcbhat_wm,WYhat_wm,WRES_wm,~,wse_wm,wtv_wm,wzv_wm] = nfast(Ywm,X,TR,tcon,pmethod,chunksize_wm); 
    
    % put back the results together.
    Wcbhat   = zeros(nvox,1);
    WYhat    = zeros(ntp,nvox);
    WRES     = zeros(ntp,nvox);
    WBLUSRES = [];
    wse      = zeros(nvox,1);
    wtv      = zeros(nvox,1);
    wzv      = zeros(nvox,1);    

    size(wse_wm), size(wse_gm) 
    
    % put them back together
    wse(Idx_wm)    = wse_wm ;    wse(Idx_gm)     = wse_gm;
    wtv(Idx_wm)    = wtv_wm ;    wtv(Idx_gm)     = wtv_gm; 
    wzv(Idx_wm)    = wzv_wm ;    wzv(Idx_gm)     = wzv_gm;
    Wcbhat(Idx_wm) = Wcbhat_wm ; Wcbhat(Idx_gm)  = Wcbhat_gm;

    WYhat(:,Idx_wm) = WYhat_wm ; WYhat(:,Idx_gm) = WYhat_gm ;
    WRES(:,Idx_wm)  = WRES_wm  ; WRES(:,Idx_gm)  = WRES_gm ;
    
end


function [cbhat,RES,stat,se,tv,zv, Wcbhat,WYhat,WRES,WBLUSRES,wse,wtv,wzv,W,V,Cy] = nfast(Y,X,TR,tcon,pmethod,chunksize)
%[W,V,Cy] = gfast(Y,X,TR)
% 
% Y:  Full time series TxV
% X:  Full design TxP 
% TR: Repetion Time [float scalar]
% 
% W:  Whitening matrix [TxT sparse matrix]
% V:  Autocorrelation Matrix [TxT matrix]
% Cy: Covariance of Y [TxT matrix]
% 
% Reimplementation of SPM FAST
% Uses p-values of overal significance for pooling
% 
% You need SPM to run this. 
% The execution time is _mainly_ dependent on the length of time series
%
% SA & TEN, Ox, 2020
%
    WBLUSRES = []; 
        
    if ~exist('pmethod','var');   pmethod = 0;       end; 
    if ~exist('chunksize','var'); chunksize = 10000; end; 

    ntp   = size(X,1); 
    nvox  = size(Y,2); 

    [cbhat,Yhat,RES,stat] = myOLS(Y,X); % to get pvalues for Fstats of overall sig.
    trRV           = stat.df; 
    
    se             = stat.se;
    tv             = stat.tval;
    zv             = stat.zval;    
    
    
    % FAST variance components for FAST
    Vi      = spm_Ce('fast',ntp,TR);
    
    if pmethod
    %%% --------------------- pooling by F-statistics
        disp(['nfast:: pooling by F-statistics.'])
        jidx  = find((stat.fp.*nvox)<0.001); % Harsh bonferroni 
        %clear stat
    else
    %%% ---------------------pooling by ACL/ACF
        disp(['nfast:: pooling by autocorrelation.'])
        [acf ,acfCI] = AC_fft(RES,ntp);
        %acl          = sum(acf(1,1:fix(ntp/4)).^2); % Anderson's suugestion of ignoring beyond ntps/4
        acf          = acf(:,1+1); %acf(1) 
        jidx         = find(abs(acf)>acfCI(2));  % only if a voxel exceeds the CI
        
       
        clear acf
    end
    
    q = numel(jidx); 
    
    if ~q
        error('nfast:: Something is wrong, there should be at least some voxels significant to the overall design'); 
    else
        disp(['nfast:: Number of pooled voxels: ' num2str(q)])
    end

    ResSS = sum(RES.^2); 

    % Whole business below is to get the (quick) effective DOF for multiple session setting. 
    % xX.xKXs      = spm_sp('Set',spm_filter(xX.K,xX.W*xX.X));
    % xX.xKXs.X    = full(xX.xKXs.X);
    % trRV  = spm_SpUtil('trRV',xX.xKXs);
    % q = spdiags(sqrt(trRV./ResSS(j)'),0,q,q) % this is HUGE! We need chunks

    %chunksize = 20000; % this is reasonable, but should be lower if low memory
    nbchunks  = ceil(q/chunksize);
    chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),q+1);

    % preallocate the variables.
    Wcbhat  = cbhat; 
    wse     = se; 
    wtv     = tv; 
    wzv     = zv; 
    WRES    = RES; 
    WYhat   = Yhat; 
    
    tic;
    for ichunk = 1:nbchunks
        chunk  = chunks(ichunk):chunks(ichunk+1)-1;
        nchunk = numel(chunk);
        disp(['nfast:: chunk ' num2str(ichunk) '/' num2str(nbchunks) ', size: ' num2str(nchunk)])
        jchunk = jidx(chunk);  
        Yc     = Y(:,jchunk);
        %sd     = 1./diag(sqrt(trRV./ResSS(jchunk)'));
        sd     = sqrt(ResSS(jchunk)/trRV);
        nYc    = Yc./sd;
        Cy     = nYc*nYc';
        Cy     = Cy./nchunk;
        
        % Call ReML to get the auto-covariance of the system
        disp('nfast:: Finding autocovariance matrix using ReML')
        V       = spm_reml(Cy,X,Vi);
        V       = V*ntp/trace(V); 
        
        % Prewhitening Matrix, W
        disp('nfast:: getting global whitening matrix.')
        W      = spm_sqrtm(spm_inv(V));
        W      = W.*(abs(W)> 1e-6);
        
        % Prewhiten the X & Y globally 
        disp('nfast:: Whiten the data and the design')
        WYc = W*Yc;
        WX  = W*X(:,2:end); %exclude the intercept while prewhitening
        
        disp('nfast:: Refit the prewhitened model.')
        WX = [ones(ntp,1) WX]; % add back the intercept
        [Wcbhatct,WYhatct,WRESct,wstat] = myOLS(WYc,WX,tcon);
        Wcbhat(jchunk)                  = Wcbhatct;
        wse(jchunk)                     = wstat.se;
        wtv(jchunk)                     = wstat.tval;
        wzv(jchunk)                     = wstat.zval;
        WRES(:,jchunk)                  = WRESct;
        WYhat(:,jchunk)                 = WYhatct;
        
    end
    toc
    disp('nfast:: put all back together.')
        
    %BLUSres this will be taking ages so we drop it for now. 
    %WBLUSRES(:,iv) = BLUSres(WY,WX,1:size(WX,2));
        
end