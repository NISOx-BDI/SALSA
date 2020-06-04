% ts_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/filtered_func_data.nii.gz';
% tcon_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design.con';
% dmat_fname='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/design_mat.txt';
% path2mask='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mask.nii.gz';
% parmat='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/FeatTest/sub-A00008326++++.feat/mc/prefiltered_func_data_mcf.par';
% %feat5/featlib.cc 
% 
% addpath('/Users/sorooshafyouni/Home/GitClone/FILM2/mis')
% addpath('/Users/sorooshafyouni/Home/matlab/spm12')
% 
% [Y,ImgStat] = CleanNIFTI_spm(ts_fname,'demean');
% Y  = Y';
% Y  = Y - mean(Y);
% T  = 900;
% TR = 0.645;
% 
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
% [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = gReMLxACF5(Y,X,TR,tcon,-2,ImgStat,path2mask,1,K);
% 
% [PSDx,PSDy]   = DrawMeSpectrum(RES,1);
% [WPSDx,WPSDy] = DrawMeSpectrum(WRES,1);
% 
% figure; hold on; grid on; 
% plot(PSDx,mean(PSDy,2))
% plot(WPSDx,mean(WPSDy,2))


function [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,WBLUSRES,wse,wtv,wzv] = gReMLxACF(Y,X,TR,tcon,tukey_m,ImgStat,path2mask,badjflag,K)
% Performs two stage prewhitening: 1) FAST 2) ACFadj
% 
% tcon should already have the intercept
% badjflag [boolean]
% K is the filter
% 
    disp(['gReMLxACF:: fit the intial naive model.'])
    % This bit naturally comes out of the gReML. But keep it here for now. 
    [cbhat,~,RES,stat] = myOLS(Y,X,tcon);
    se                 = stat.se;
    tv                 = stat.tval;
    zv                 = stat.zval;

    disp(['gReMLxACF:: fit gloabl FAST.'])
    
    [WY,WX]            = gReML(Y,X,TR,tcon,0); 
    
    disp(['gReMLxACF:: fit voxel-wise prewhitening.'])
    
    [~,~,~,~,~,~,Wcbhat,WYhat,WRES,WBLUSRES,wse,wtv,wzv] = feat5(WY,WX,tcon,tukey_m,ImgStat,path2mask,badjflag,0,K);
  
    
end


function [WY,WX,RES] = gReML(Y,X,TR,tcon,pmethod)
%[W,V,Cy] = gReML(Y,X,TR)
% This is a ligher version of gfast.m
% 
% Y:  Full time series TxV
% X:  Full design TxP 
% TR: Repetion Time [float scalar]
% 
% WY:  Whitened Y
% WY:  Whitened X 
% 
% Reimplementation of SPM FAST
% Uses p-values of overal significance for pooling
% 
% You need SPM to run this. 
% The execution time is _mainly_ dependent on the length of time series
%
% SA & TEN, Ox, 2020
%

    if nargin<5; pmethod = 0; end; 

    ntp   = size(X,1); 
    nvox  = size(Y,2); 

    [~,~,RES,stat] = myOLS(Y,X,tcon); % to get pvalues for Fstats of overall sig.
    trRV           = stat.df; 
        
    if pmethod
    %%% --------------------- pooling by F-statistics
        disp(['gReML:: pooling by F-statistics.'])
        jidx  = find((stat.fp.*nvox)<0.001); % Harsh bonferroni 
        clear stat
    else
    %%% ---------------------pooling by ACL/ACF
        disp(['gReML:: pooling by autocorrelation.'])
        [acf ,acfCI] = AC_fft(RES,ntp);
        %acl          = sum(acf(1,1:fix(ntp/4)).^2); % Anderson's suugestion of ignoring beyond ntps/4
        acf          = acf(:,1+1); %acf(1) 
        jidx         = find(abs(acf)>acfCI(2));  % only if a voxel exceeds the CI
        
        clear acf
    end
    
    q = numel(jidx); 
    
    if ~q
        error('gReML:: Something is wrong, there should be at least some voxels significant to the overall design'); 
    else
        disp(['gReML:: Number of pooled voxels: ' num2str(q)])
    end

    ResSS = sum(RES.^2); 

    % Using chunks is kinda redundent in matlab because of broadcasting. 
    % But, it is better to keep it this way in case we wanted to implement
    % this later in C++ OR python. 
    % 
    chunksize = 2000; % this is reasonable, but should be lower if low memory
    nbchunks  = ceil(q/chunksize);
    chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),q+1);

    Cy = 0; 
    for ichunk = 1:nbchunks 
        disp(['gReML:: chunk ' num2str(ichunk) '/' num2str(nbchunks)])
        chunk  = chunks(ichunk):chunks(ichunk+1)-1;
        jchunk = jidx(chunk);  
        %sd     = 1./diag(sqrt(trRV./ResSS(jchunk)'));
        sd     = sqrt(ResSS(jchunk)/trRV);
        Yc     = Y(:,jchunk)./sd;
        Cy     = Cy + Yc*Yc';
    end
    Cy = Cy/q; %Average across the pool 

    clear Yc v 

    % Error variance components for FAST
    J        = 0; % Keep the basis limited to AR(1)
    Vi       = CovFast(ntp,TR,J);
    %Vi      = spm_Ce('fast',ntp,TR);

    % Call ReML to get the auto-covariance of the system
    disp('gReML:: Finding autocovariance matrix using ReML')
    V       = spm_reml(Cy,X,Vi);
    V       = V*ntp/trace(V); 

    % Prewhitening Matrix, W
    disp('gReML:: getting global whitening matrix.')
    W      = spm_sqrtm(spm_inv(V));
    W      = W.*(abs(W)> 1e-6);
    
    % Prewhiten the X & Y globally 
    disp('gReML:: Prewhiten the data and the design.')
    WY = W*Y;
    WX = W*X(:,2:end); %exclude the intercept while prewhitening
    WX = [ones(ntp,1), WX];
end



function C = CovFast(ntp,TR,J)
% A version of spm_Ce.m with 'fast' option but with control over the basis
% ntp: # data points [integer]
% TR : repetion time [float]
% J  : # of basis beyond AR(1) [integer]
% NB! The function can't deal with multiple sessions.
% SA, Ox, 2020
if nargin<3; J = 2; end; % if not specified, do what spm_Ce.m does
    C  = {};
    T     = (0:(ntp - 1))*TR;
    d     = 2.^(floor(log2(TR/4)):log2(64));
    for i = 1:min(6,length(d))
        for j = 0:J
            QQ = toeplitz((T.^j).*exp(-T/d(i)));
            [x,y,q] = find(QQ);
            C{end + 1} = sparse(x,y,q,ntp,ntp);
        end
    end

end

