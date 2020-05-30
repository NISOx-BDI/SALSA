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
% W = gfast5(Y,X,0.645); 



function [W,V,Cy] = gfast(Y,X,TR)
% Reimplementation of SPM FAST
% Uses p-values of overal significance for pooling
% 
% You need SPM to run this. 
% The execution time is _mainly_ dependent of the time series length not the number
% of voxels. 
%
% SA, Ox, 2020
%

ntp   = size(X,1); 
nvox  = size(Y,2); 

[~,~,RES,stat] = myOLS(Y,X); % to get pvalues for Fstats of overall sig.
jidx  = find((stat.fp.*size(Y,2))<0.001); % Harsh bonferroni 
q     = numel(jidx); 

if ~q
    error('gfast:: Something is wrong, there should be at least some voxels significant to the overal design'); 
else
    disp(['gfast:: Number of pooled voxels: ' num2str(q)])
end

ResSS = sum(RES.^2); 

% Whole business below is to get the (quick) effective DOF for multiple session setting. 
% xX.xKXs      = spm_sp('Set',spm_filter(xX.K,xX.W*xX.X));
% xX.xKXs.X    = full(xX.xKXs.X);
% trRV  = spm_SpUtil('trRV',xX.xKXs);
% q = spdiags(sqrt(trRV./ResSS(j)'),0,q,q) % this is HUGE! We need chunks

trRV      = stat.df; clear stat RES
chunksize = 2000; % this is reasonable, but should be lower if low memory
nbchunks  = ceil(q/chunksize);
chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),q+1);


Cy = 0; 
for ichunk = 1:nbchunks
    disp(['gfast:: chunk ' num2str(ichunk) '/' num2str(nbchunks)])
    chunk  = chunks(ichunk):chunks(ichunk+1)-1;
    jchunk = jidx(chunk);  
    v      = diag(sqrt(trRV./ResSS(jchunk)'));
    
    Yc     = Y(:,jchunk)*v;
    Cy     = Cy + Yc*Yc';
end
Cy = Cy/q;

clear Y Yc v 

% FAST variance components for FAST
Vi      = spm_Ce('fast',ntp,TR);

% Call ReML to get the auto-covariance of the system
disp('gfast:: Finding autocovariance matrix using ReML')
V       = spm_reml(Cy,X,Vi);
V       = V*ntp/trace(V); % make sure the covariance matrix is not non-degenrate

% Prewhitening Matrix, W
disp('gfast:: getting global whitening matrix.')
W      = spm_sqrtm(spm_inv(V));
W      = W.*(abs(W)> 1e-6);

end