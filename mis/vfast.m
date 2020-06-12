% clear
% 
% clear
% 
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
% T=900;
% TR = 0.645; 
% 
% disp('MC params.')
% MCp      = load(parmat); 
% MCp      = GenMotionParam(MCp,24); 
% X        = [load(dmat_fname) MCp];
% 
% X        = X-mean(X); 
% 
% disp('hpf')
% K = hp_fsl(size(Y,1),100,TR);    
% X     = K*X;    % high pass filter the design
% Y     = K*Y;  % high pass filter the data
% 
% X = [ones(T,1) X];
% tcon     = zeros(1,size(X,2));
% tcon(2)  = 1;
% 
% % [cbhat_gm,RES_gm,stat_gm,se_gm,tv_gm,zv_gm,Wcbhat_gm,WYhat_gm,WRES_gm,wse_gm,wtv_gm,wzv_gm] = nfast5(Ygm,X,0.645,tcon); 
% % [PSDx_gm,PSDy_gm]   = DrawMeSpectrum(RES_gm,1);
% % [WPSDx_gm,WPSDy_gm] = DrawMeSpectrum(WRES_gm,1);
% % [cbhat_wm,RES_wm,stat_wm,se_wm,tv_wm,zv_wm,Wcbhat_wm,WYhat_wm,WRES_wm,wse_wm,wtv_wm,wzv_wm] = nfast5(Ywm,X,0.645,tcon); 
% % [PSDx_wm,PSDy_wm]   = DrawMeSpectrum(RES_wm,1);
% % [WPSDx_wm,WPSDy_wm] = DrawMeSpectrum(WRES_wm,1);
% % figure; hold on; grid on; 
% % plot(PSDx_wm,mean(PSDy_wm,2))
% % plot(WPSDx_wm,mean(WPSDy_wm,2))
% 
% ARO = 20; 
% [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = vfast5(Y,X,TR,tcon,ARO);
% 
% [PSDx,PSDy]   = DrawMeSpectrum(RES,1);
% [WPSDx,WPSDy] = DrawMeSpectrum(WRES,1);
% 
% figure; hold on; grid on;
% plot(PSDx,mean(PSDy,2))
% plot(WPSDx,mean(WPSDy,2))

function [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = vfast(Y,X,TR,tcon,ARO)


[ntp,nvox]      = size(Y);

% Only once for a brain-----------------------
disp(['vfast:: getting residuals.'])
pinvX           = pinv(X); 
R               = eye(ntp)-X*pinvX; % residual forming matrix 
RES             = R*Y; 

[cbhat,~,~,stat]  = myOLS(Y,X,tcon); % just for the evaluation, remove later.
se                = stat.se;
tv                = stat.tval;
zv                = stat.zval;
% Get the autocovariances
disp(['vfast:: getting autocovariances.'])
[~,~,autocovs]  = AC_fft(RES,ntp);
a               = autocovs(:,1:ARO+1)';

% matrix of basis
disp(['vfast:: getting matrix of AR basis.'])
Bfull           = ACFbasis(ntp-1,TR); 
MaxLag          = ARO+1;  % I don't like calling this ARO, but, whatever we call it
                          % it is the key free parameter
B               = Bfull(1:MaxLag,:);

% matrix of M
disp(['vfast:: Autocovariance adjustment matrix.'])
M               = BiasAdjMat(R,ntp,ARO);

disp(['vfast:: getting FAST basis coefficients.'])
g              = (M*B)\a;
g              = g./g(1,:); % normalise the FAST basis, assuming the same is valid here as it is in Appendix A

% Pervoxel -----------------------------------
% 
% pre-allocate memory
WYhat  = zeros(ntp,nvox); WRES = WYhat;
Wcbhat = zeros(1,nvox); wse = Wcbhat; wtv = wse; wzv = wtv;
for iv = 1:nvox
    if ~mod(iv,5000); disp(['on voxel: ' num2str(iv)]); end; 
    
    % Worsely's quick W 
    [Ainvt,posdef]  = chol(toep(g(:,iv))); 
    p1              = size(Ainvt,1); % this is basically posdef - 1, but we'll keep it as Keith Worsely's code. 
    A               = inv(Ainvt'); 
    
    W               = toep([A(p1,p1:-1:1) zeros(1,ntp-p1)],zeros(1,ntp)); 
    W(1:p1,1:p1)    = A;    
    W               = W.*(abs(W)> 1e-6); % not really important if I use Worsely's code. 
    
    % Whiten X & Y
    WYv        = W*Y(:,iv);
    WX         = W*X;
    
    % Refit the model
    [Wcbhat(iv),WYhat(:,iv),WRES(:,iv),wstat] = myOLS(WYv,WX,tcon); 
    wse(iv)  = wstat.se;
    wtv(iv)  = wstat.tval;
    wzv(iv)  = wstat.zval;
end

 function B = ACFbasis(MaxLag,TR)
 % matrix of Mxlag x basis
 % TEN, 2020
    t = [(0:MaxLag)*TR]';
    d = 2.^(floor(log2(TR/4)):log2(64));
    d(7:end) = [];
    j = [0 1];%[0 1]; % or	[0 1 2], but I don't think j=2 buys us much
    J = repmat(j,length(t),length(d));
    T = repmat(t,1,length(j)*length(d));
    D = repmat(repelem(d,length(j)),length(t),1);
    B = T.^J.*exp(-T./D);
    B = B./std(B);
 end

 function M = BiasAdjMat(R,ntp,ARO)
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
    M   = zeros(ARO+1);
    for l = 1:(ARO+1)
        Dl = diag(ones(1,ntp-l+1),l-1); % upper triangle
        for j = 1:(ARO+1)
           DjDjt    = (diag(ones(1,ntp-j+1),j-1)+diag(ones(1,ntp-j+1),-j+1))/(1+(j==1));
           M(l,j)   = trace(R'*Dl*R*DjDjt);
        end
    end
% ------------------------------------------------------------------------

 end
 
 function t = toep(c,r)
    % quick toeplitz to get rid of the checks and warnings
    %Toeplitz matrix.
    if nargin < 2
        c(1) = conj(c(1));                    % set up for Hermitian Toeplitz
        r = c; 
        c = conj(c); 
    end
    r = r(:);                               % force column structure
    c = c(:);
    p = length(r);
    m = length(c);
    x = [r(p:-1:2, 1) ; c];                 % build vector of user data
    ij = (0:m-1)' + (p:-1:1);               % Toeplitz subscripts
    t = x(ij);                              % actual data
    if isrow(ij)                            % preserve shape for a single row
        t = t.';
    end
 end

end