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
% Y = Y';
% Y = Y - mean(Y);
% T=900;
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
% [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,wse,wtv,wzv] = arw5(Y,X,tcon,20,ImgStat,[],[]);
% 
% [PSDx,PSDy]   = DrawMeSpectrum(RES,1);
% [WPSDx,WPSDy] = DrawMeSpectrum(WRES,1);
% 
% figure; hold on; grid on; 
% plot(PSDx,mean(PSDy,2))
% plot(WPSDx,mean(WPSDy,2))


function [cbhat,RES,stat,se,tv,zv,Wcbhat,WYhat,WRES,WBLUSRES,wse,wtv,wzv] = arw5(Y,X,tcon,ARO,ImgStat,path2mask,K)
% Y    : TxV
% X    : TxEV. X should have the detreding basis + motion parameters
% tcon : 1xEV
%
%
% SA, Ox, 2020

    WBLUSRES = []; 
    
    warning('off','all')

    disp('::arw::')

    ntp                   = size(Y,1); 
    nvox                  = size(Y,2); 
    [cbhat,~,~,stat] = myOLS(Y,X,tcon); % if cbhat is not needed; R=I-(XX+); RES=RY

    se                   = stat.se;
    tv                   = stat.tval;
    zv                   = stat.zval;
    
    pinvX               = pinv(X); 
    ResidFormingMat     = eye(ntp)-X*pinvX; % residual forming matrix 
    RES                 = ResidFormingMat*Y;

    if ~exist('K','var')  || isempty(K)
        K = eye(ntp);   
    end;     
    
    % find the acf of the residuals
    [~,~,dRESacov]      = AC_fft(RES,ntp); % Autocovariance; VxT
    dRESacov            = dRESacov';
    dRESacov            = dRESacov(1:ARO+1,:); %cut-off here so smoothing gets much faster
    
    % smooth ACF
    if ~isempty(path2mask)
        acfFWHMl               = 5; 
        dRESacov            = ApplyFSLSmoothing(dRESacov',acfFWHMl,ImgStat,path2mask)';
        disp(['arw:: Estimate ACF, smooth on ' num2str(acfFWHMl) 'mm and taper on ' num2str(ARO) ' lag.'])
    end
    % bias adjusting matrix + filter K 
    BiasAdj             = ACFBiasAdjMat(ResidFormingMat*K,ntp,ARO); 

    
    % refit the model to pre-whitened data
    disp('arw:: Refit the prewhitened model.')
    Wcbhat   = zeros(1,nvox);
    WYhat    = zeros(ntp,nvox); 
    WRES     = zeros(ntp,nvox);
    %WBLUSRES = zeros(ntp,nvox);
    wse      = zeros(nvox,1); 
    wtv      = zeros(nvox,1);
    wzv      = zeros(nvox,1);
    for iv = 1:nvox
        if ~mod(iv,5000); disp(['on voxel: ' num2str(iv)]); end; 
        sqrtmVhalf = establish_prewhiten_matrix(dRESacov(:,iv),ntp,BiasAdj);
        WYv        = sqrtmVhalf*Y(:,iv);
        WX         = sqrtmVhalf*X;   
        [Wcbhat(iv),WYhat(:,iv),WRES(:,iv),wstat] = myOLS(WYv,WX,tcon);
        
        wse(iv)  = wstat.se;
        wtv(iv)  = wstat.tval;
        wzv(iv)  = wstat.zval;
        %WY(:,iv)   = WYv;
        
        % BLUSres -- this will take ages. 
        % WBLUSRES(:,iv)                            = BLUSres(WYv,WX,1:size(WX,2));        
    end

end

function [sqrtmVhalf,posdef] = establish_prewhiten_matrix(autocov,ntp,BiasAdj)
% dRESacov : autocovariance 
% Mord : Order required 
% BiasAdj: bias adjustment matrix 
%
% SA, Ox, 2020
%
     
    %dRESac_adj              = (BiasAdj*autocov(1:AROrd+1));
    dRESac_adj              = BiasAdj*autocov;
    dRESac_adj              = dRESac_adj./dRESac_adj(1); % make a auto*correlation*
    
    [Ainvt,posdef]          = chol(toep(dRESac_adj)); 
    p1                      = size(Ainvt,1); % this is basically posdef - 1, but we'll keep it as Keith Worsely's code. 
    A                       = inv(Ainvt'); 
    
    sqrtmVhalf              = toep([A(p1,p1:-1:1) zeros(1,ntp-p1)],zeros(1,ntp)); 
    sqrtmVhalf(1:p1,1:p1)   = A;
    
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