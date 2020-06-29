% ts_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/filtered_func_data.nii.gz';
% tcon_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design.con';
% dmat_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design_mat.txt';
% path2mask='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mask.nii.gz';
% parmat='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mc/prefiltered_func_data_mcf.par';
% %feat5/featlib.cc 
% 
% addpath('/Users/sorooshafyouni/Home/GitClone/FILM2/mis')
% addpath('/Users/sorooshafyouni/Home/GitClone/FILM2/utils/Trend')
% addpath('/Users/sorooshafyouni/Home/matlab/spm12')
% 
% [Y,ImgStat] = CleanNIFTI_spm(ts_fname,'demean');
% Y = Y';
% Y = Y - mean(Y);
% T=900;
% 
% disp('MC params.')
% MCp      = load(parmat); 
% MCp      = GenMotionParam(MCp,24); 
% X        = [load(dmat_fname) MCp];
% 
% X        = X-mean(X); 
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
% 
% [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gfast5(Y,X,0.645,tcon); 
% 
% [PSDx,PSDy]   = DrawMeSpectrum(RES,1);
% [WPSDx,WPSDy] = DrawMeSpectrum(WRES,1);
% 
% figure; hold on; grid on; 
% plot(PSDx,mean(PSDy,2))
% plot(WPSDx,mean(WPSDy,2))

function [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,WBLUSRES,wse,wtv,wzv,W,V,Cy] = gfast(Y,X,TR,tcon,pmethod)
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
    
    if nargin<5; pmethod = 0; end; 

    ntp   = size(X,1); 
    nvox  = size(Y,2); 

    [cbhat,~,RES,stat] = myOLS(Y,X); % to get pvalues for Fstats of overall sig.
    trRV           = stat.df; 
    
    se             = stat.se;
    tv             = stat.tval;
    zv             = stat.zval;    
    
    if pmethod
    %%% --------------------- pooling by F-statistics
        disp(['gfast:: pooling by F-statistics.'])
        jidx  = find((stat.fp.*nvox)<0.001); % Harsh bonferroni 
        %clear stat
    else
    %%% ---------------------pooling by ACL/ACF
        disp(['gfast:: pooling by autocorrelation.'])
        [acf ,acfCI] = AC_fft(RES,ntp);
        %acl          = sum(acf(1,1:fix(ntp/4)).^2); % Anderson's suugestion of ignoring beyond ntps/4
        acf          = acf(:,1+1); %acf(1) 
        jidx         = find(abs(acf)>acfCI(2));  % only if a voxel exceeds the CI
        
       
        clear acf
    end
    
    q = numel(jidx); 
    
    if ~q
        error('gfast:: Something is wrong, there should be at least some voxels significant to the overall design'); 
    else
        disp(['gfast:: Number of pooled voxels: ' num2str(q)])
    end

    ResSS = sum(RES.^2); 

    % Whole business below is to get the (quick) effective DOF for multiple session setting. 
    % xX.xKXs      = spm_sp('Set',spm_filter(xX.K,xX.W*xX.X));
    % xX.xKXs.X    = full(xX.xKXs.X);
    % trRV  = spm_SpUtil('trRV',xX.xKXs);
    % q = spdiags(sqrt(trRV./ResSS(j)'),0,q,q) % this is HUGE! We need chunks

    chunksize = 2000; % this is reasonable, but should be lower if low memory
    nbchunks  = ceil(q/chunksize);
    chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),q+1);


    Cy = 0; 
    for ichunk = 1:nbchunks
        disp(['gfast:: chunk ' num2str(ichunk) '/' num2str(nbchunks)])
        chunk  = chunks(ichunk):chunks(ichunk+1)-1;
        jchunk = jidx(chunk);  
        %sd     = 1./diag(sqrt(trRV./ResSS(jchunk)'));
        sd     = sqrt(ResSS(jchunk)/trRV);
        Yc     = Y(:,jchunk)./sd;
        Cy     = Cy + Yc*Yc';
    end
    Cy = Cy/q; %Average across the pool 

    clear Yc v 

    % FAST variance components for FAST
    Vi      = spm_Ce('fast',ntp,TR);

    tic; 
    % Call ReML to get the auto-covariance of the system
    disp('gfast:: Finding autocovariance matrix using ReML')
    V       = spm_reml(Cy,X,Vi);
    V       = V*ntp/trace(V); 
    
    % Prewhitening Matrix, W
    disp('gfast:: getting global whitening matrix.')
    W      = spm_sqrtm(spm_inv(V));
    W      = W.*(abs(W)> 1e-6);
    
    toc
    
    % Prewhiten the X & Y globally 
    disp('gfast:: Whiten the data and the design')
    WY = W*Y;
    WX = W*X(:,2:end); %exclude the intercept while prewhitening
    
    disp('gfast:: Refit the prewhitened model.')
    WX = [ones(ntp,1) WX]; % add back the intercept
    [Wcbhat,WYhat,WRES,wstat] = myOLS(WY,WX,tcon);
    wse  = wstat.se;
    wtv  = wstat.tval;
    wzv  = wstat.zval;
    
    %BLUSres this will be taking ages so we drop it for now. 
    %WBLUSRES(:,iv) = BLUSres(WY,WX,1:size(WX,2));
        
end