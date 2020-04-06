%clear

warning('off','all')

pwdmethod = 'ARMAHR'; %ACF AR-YW AR-W ARMAHR
Mord      = 5; 

SubID     = 'A00029304';

SesID     = 'DS2';
TR        = 0.645;

PATH2AUX='/Users/sorooshafyouni/Home/GitClone/FILM2';
addpath([PATH2AUX '/utils/Trend'])
addpath('/Users/sorooshafyouni/Home/matlab/spm12')
addpath([PATH2AUX '/utils/AR_YW'])
addpath([PATH2AUX '/utils/ARMA_HR'])
addpath([PATH2AUX '/mis'])

%Raw Images (MMP feat output)
Path2ImgRaw = [COHORTDIR '/R_mpp'];
Path2ImgDir = [Path2ImgRaw '/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-' num2str(TR*1000) '_bold_mpp'];
Path2Img    = [Path2ImgDir '/prefiltered_func_data_bet.nii'];
Path2MC     = [Path2ImgDir '/prefiltered_func_data_mcf.par'];

% Directory 2 save the results
Path2ImgResults=[COHORTDIR '/R.PW/' pwdmethod '_AR-' num2str(Mord) '/RNullfMRI_' SubID '_' SesID];

SaveImagesFlag      = 1; 
DoDetrendingPrior   = 0; 
MParamNum           = 6; 
NumTmpTrend         = 3;

%%% Read The Data %%%%%%%%%%%%%%%%%%%%%%%%
[Y,InputImgStat]=CleanNIFTI_spm(Path2Img,'demean');
T = InputImgStat.CleanedDim(2);
TR = InputImgStat.voxelsize(4);
Vorig = InputImgStat.CleanedDim(1);
V = Vorig;

%%% Generate a Design Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EDtype = 'boxcar'; 
BCl = 20;
EDX = GenerateED(BCl,T,TR); 
EDX = EDX - mean(EDX); 

%%% DETREND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DoDetrendingPrior
    pp = 1+fix(TR.*T/150);
    dY = multpolyfit(repmat(1:T,Vorig,1),Y,T,pp);
    dY = dY - mean(dY,2); 
else
    dY = Y; 
end

%%% SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Order of AR which will be used to simulated stuff
% ARporder = 20;
% ARParam  = AR_YW(dYorig,ARporder);
% nRlz = 1000; 
% model = arima('Constant',0,'AR',ARParam,'Variance',1);
% dY = simulate(model,T,'NumPaths',nRlz)'; 
% dY = dY - mean(dY,2); 

%%% DESIGN MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal trends
TempTrend = []; MCp = [];

if NumTmpTrend>0
    [TempTrend,NUMSPLINES]   = GenerateTemporalTrends(T,TR,NumTmpTrend); % DC + Poly trends + Spline trends 
    TempTrend   = TempTrend(:,2:end); % we add a column of one later.
end
disp(['Number of temporal trends: ' num2str(NumTmpTrend)])

% Motion parameters 
if MParamNum     == 6
    MCp = load(Path2MC);
elseif MParamNum == 12
    MCp = load(Path2MC);
    MCp = [MCp,MCp.^2]; % 12 parameter motion 
elseif MParamNum == 24
    % SHOULD BE DONE LATER
end
disp(['Number of motion parameter: ' num2str(MParamNum)])
    
%
X           = [EDX,TempTrend,MCp];
X           = X - mean(X); % demean everything 
ED_IDX      = 2; % where will be the expermintal design in the final X?
%%% RESIDUALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pinvX           = pinv(X); 
ResidFormingMat = eye(T)-X*pinvX; % residual forming matrix 

residY = zeros(V,T);
for vi = 1:V    
    if ~mod(vi,10000); disp(['Residuals ::: on voxel ' num2str(vi)]); end;
    residY(vi,:) = ResidFormingMat*dY(vi,:)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BIAS REDUCTION OF AUTOREGRESSIVE MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

YWflag = 0; WrosleyFlag = 0; ACFflag = 0; ARMAHRflag = 0; MPparamNum = 0; 
if strcmpi(pwdmethod,'AR-W') %Worsely %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WrosleyFlag = 1; 
    warning('off','MATLAB:toeplitz:DiagonalConflict')
    M_biasred   = zeros(Mord+1);
    for i=1:(Mord+1)
        Di                  = (diag(ones(1,T-i+1),i-1)+diag(ones(1,T-i+1),-i+1))/(1+(i==1));
        for j=1:(Mord+1)
           Dj               = (diag(ones(1,T-j+1),j-1)+diag(ones(1,T-j+1),-j+1))/(1+(j==1));
           M_biasred(i,j)   = trace( ResidFormingMat*Di*ResidFormingMat*Dj )/(1+(i>1));
        end
    end
    invM_biasred = inv(M_biasred);
elseif strcmpi(pwdmethod,'AR-YW') % Yule-Walker %%%%%%%%%%%%%%%%%%%%%%%%%%%
    YWflag  = 1; 
    invM_biasred = eye(Mord+1);
elseif strcmpi(pwdmethod,'ACF') % ACF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ACFflag = 1;
elseif strcmpi(pwdmethod,'ARMAHR') % ARMAHR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ARMAHRflag = 1; 
    MPparamNum          = 1;
    ARParamARMA         = 50;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PREWHITEN THE RESIDULAS & ESTIMATE BIAS AND CPS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% ACFs %%%%%%%%%%%%%%%
disp(['Calculate the autocorrelation coefficients.'])
[~,~,dRESacov]  = AC_fft(residY,T); % Autocovariance 
dRESacorr       = dRESacov./sum(abs(residY).^2,2); % Autocorrelation
ACL             = sum(dRESacorr.^2,2); % Autocorrelation Length

for vi = 1:V
    
    if YWflag % Yule-Walker %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp(['AR-YW ::: on voxel ' num2str(vi)]); end; 
        % Prewhiten the residuals using Yule-Walker estimates
        %YWARparam_tmp   = AR_YW(dYorig,Mord);
        %YWARparam(vi,:) = YWARparam_tmp;
       
        %AR_YW -------------------------------------------
        ac_tmp          = dRESacorr(vi,:);    
        R_tmp           = toeplitz(ac_tmp(1:Mord));
        r_tmp           = ac_tmp(2:Mord+1);
        YWARparam_tmp   = pinv(R_tmp)*r_tmp';
        % ------------------------------------------------
        
        ACMat           = full(spm_Q(YWARparam_tmp,T));
        
        invACMat        = inv(ACMat); % pinv is damn slow!
        sqrtmVhalf      = chol(invACMat);
        
    elseif ACFflag % ACF - Tukey tapered %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp(['ACF-TUKEY ::: on voxel ' num2str(vi)]); end; 
        acfdRES     = dRESacorr(vi,:); 
        %Mord        = 2*round(sqrt(T));
        acfdRES_tt  = [1 TukeyTaperMe(acfdRES(2:end),T-1,Mord)];
        ACMat       = toeplitz(acfdRES_tt);  
        
        invACMat   = inv(ACMat); % pinv is damn slow!
        sqrtmVhalf = chol(invACMat);
        
    elseif WrosleyFlag % Worsely %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp(['AR-W ::: on voxel ' num2str(vi)]); end; 
        dRESac_adj      = (invM_biasred*dRESacov(vi,1:Mord+1)');
        dRESac_adj      = dRESac_adj./dRESac_adj(1); % make a auto*correlation*
        
        [Ainvt,posdef]      = chol(toeplitz(dRESac_adj)); 
        p1                  = size(Ainvt,1); 
        A                   = inv(Ainvt'); 
        sqrtmVhalf          = toeplitz([A(p1,p1:-1:1) zeros(1,T-p1)],zeros(1,T)); 
        sqrtmVhalf(1:p1,1:p1) = A;
        
        
    elseif ARMAHRflag % ARMA HR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp(['ARMA-HR ::: on voxel ' num2str(vi)]); end; 
        
        %AR_YW -------------------------------------------
        ac_tmp              = dRESacorr(vi,:);    
        R_tmp               = toeplitz(ac_tmp(1:ARParamARMA));
        r_tmp               = ac_tmp(2:ARParamARMA+1);
        YWARparam_tmp       = pinv(R_tmp)*r_tmp';
        % ------------------------------------------------
        
        [arParam,maParam]   = ARMA_HR_ACF(residY(vi,:)',YWARparam_tmp',T,Mord,MPparamNum);
        %[arParam,maParam]   = ARMA_HR(residY(vi,:)',Mord,MPparamNum);
        ACMat               = ARMACovMat([arParam,maParam],T,Mord,MPparamNum);

        invACMat   = inv(ACMat); % pinv is damn slow!
        sqrtmVhalf = chol(invACMat);        
    end
        
    % Make the X & Y whitened %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ystar_YW = sqrtmVhalf*dY(vi,:)';
    Xstar_YW = sqrtmVhalf*X;
    
    % Fit a model to the prewhitened system  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    Xstar_YW                                         = [ones(T,1), Xstar_YW]; % add intercept
    [Bhat_PW_S_tmp,~,dpwRES_tmp,Stat_PW_SE_T_tmp]    = myOLS(Ystar_YW,Xstar_YW);
    Bhat_PW(vi)    = Bhat_PW_S_tmp;         
    SE_PW(vi)      = Stat_PW_SE_T_tmp.se;
    tVALUE_PW(vi)  = Stat_PW_SE_T_tmp.tval;
    
    % Whithout prewhitening of error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X0                                                  = [ones(T,1),X];
    [Bhat_Naive_S_tmp,~,resNaive,Stat_Naive_SE_tmp]   = myOLS(dY(vi,:),X0);
    Bhat_Naive(vi)    = Bhat_Naive_S_tmp;
    SE_Naive(vi)      = Stat_Naive_SE_tmp.se;
    tVALUE_Naive(vi)  = Stat_Naive_SE_tmp.tval;
    
    [~,CPSstat_PW(vi),CPZ_PW(vi)]       = CPSUnivar(dpwRES_tmp,X0);
    [~,CPSstat_Naive(vi),CPZ_Naive(vi)] = CPSUnivar(resNaive,X0);
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE THE RESULTS AS AN IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SaveImagesFlag
    % 3D IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VariableList = {'Bhat_Naive','SE_Naive','tVALUE_Naive',...
        'Bhat_PW','SE_PW','tVALUE_PW',...
        'CPSstat_PW','CPZ_PW',...
        'CPZ_Naive','CPSstat_Naive'};
    OutputImgStat            = InputImgStat.spmV(1);
    OutputImgStat.Removables = InputImgStat.Removables;

    for vname = VariableList

        tmpvar                   = eval(vname{1});
        OutputImgStat.fname      = [Path2ImgResults '/ED' EDtype '_' num2str(BCl) '_' pwdmethod '_AR' num2str(Mord) '_MA' num2str(MPparamNum)  '_' vname{1} '.nii'];

        CleanNIFTI_spm(tmpvar,'ImgInfo',OutputImgStat);

    end
end

% 4D IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% Beta_PW_SE_S_STD    = std(Beta_PW_SE_S); 
% Beta_PW_SE_T_MEAN   = mean(Beta_PW_SE_T);
% 
% Bhat_NAIVE_S_STD      = std(Bhat_Naive_S); 
% Beta_NAIVE_SE_T_MEAN  = mean(Beta_Naive_SE_T); 

% save(['R/NullSim_Univar_Biases_' pwdmethod '_AR' num2str(Mord) '.mat'],...
%     'Beta_PW_SE_S','rngid','Beta_PW_SE_T','Bhat_Naive_S','Beta_Naive_SE_T',...
%     'CPSstat_Naive','CPZ_Naive','CPSstat_PW','CPZ_PW',...
%     'dY','Y','T','TR')

% BetaBias = @(Truth_MC,Estimated_TH) ((Estimated_TH-Truth_MC)./Truth_MC).*100;
% fh = figure; 
% hold on; grid on; box on; 
% title('[THEORETICAL - (TRUE) SIMULATION]/ (TRUE) SIMULATION x 100')
% bh = bar([BetaBias(Bhat_NAIVE_S_STD,Beta_NAIVE_SE_T_MEAN),BetaBias(Beta_PW_SE_S_STD,Beta_PW_SE_T_MEAN)]);
% fh.Children.XTick=1:2;
% fh.Children.XTickLabel = {'Naive','Prewhitened'};
% ylabel('Bias in SE of $\hat\beta$','Interpreter','latex')


%Sanity check
% for i = 1:1000
%     Y=randn(1,900); 
%     Bhat_tmp = X\Y'; sanitycheck: [b,c,d] = glmfit(Xdpw,YYYdpw);
%     Yhat = X*Bhat_tmp;
%     res  = Y'-Yhat;
%     
%     Bhat(i) = Bhat_tmp;
% 
%     Bhat_SE(i) = sqrt(((res'*res))./(X'*X)/(T-2));
% end

% if WrosleyFlag
%     M_biasred=zeros(ARord+1);
%     for i=1:(ARord+1)
%         Di=(diag(ones(1,T-i+1),i-1)+diag(ones(1,T-i+1),-i+1))/(1+(i==1));
%         for j=1:(ARord+1)
%            Dj=(diag(ones(1,T-j+1),j-1)+diag(ones(1,T-j+1),-j+1))/(1+(j==1));
%            M_biasred(i,j)=trace(ResidFormingMat*Di*ResidFormingMat*Dj)/(1+(i>1));
%         end
%     end
%     invM_biasred = inv(M_biasred);
% else
%     invM_biasred = eye(ARord+1);
% end
% 
% dpwRES      = zeros(V,T); 

%%% PCA on the RAW Data %%%%%%%%%%%%%%%%%%
% the data has already been dmeaned inside CleanNIFTI function. 
% [EigCovMat,PCs,LTfac,~,VarExp,mu]=pca(Y');
% 
% fh_pca = figure('position',[50,500,1000,600]); 
% hold on; 
% for icp=1:10
%     subplot(5,2,icp); hold on; grid on; box on
%     title(['PC ' num2str(icp) ' , VE: ' num2str(VarExp(icp))])
%     plot(PCs(:,icp)); xlim([0 T])
% end
% set(fh_pca,'color','w')

%50,9000,4000,35000
%idxY = 35000;
% dt_examples_fh = figure('position',[50,500,600,600]);
% ii = 1; 
% for idxY = [50,9000,4000,35000]
%     subplot(2,2,ii)
%     hold on; box on; grid on; 
%     title('Detrending')
%     plot(Y(idxY,:));
%     plot(dY(idxY,:));
%     xlim([0 T])
%     legend({'Raw BOLD','Deterended BOLD'})
%     ylabel('a.u.'); xlabel('Scan')
%     ii = ii + 1; 
% end
% set(dt_examples_fh,'color','w')
% export_fig(dt_examples_fh,'RFig/dY.png')

%%% PCA on the deterended Data %%%%%%%%%%%
% [dEigCovMat,dPCs,dLTfac,~,dVarExp,dmu]=pca(dY');
% 
% fh_dpca = figure('position',[50,500,1000,600]); 
% hold on; 
% for icp=1:10
%     subplot(5,2,icp); hold on; grid on; box on
%     title(['PC ' num2str(icp) ' , VE: ' num2str(dVarExp(icp))])
%     
%     plot(dPCs(:,icp)); xlim([0 T])
% end
% set(fh_dpca,'color','w')

%%% ESTIMATE AR %%%%%%%%%%%%%%%%%%%%%%%%%%
% arp = 20; 
% dARp = AR_YW_voxel(dY,T,arp);

% [ARdEigCovMat,ARdPCs,ARdLTfac,~,ARdVarExp,ARdmu]=pca(dARp'-mean(dARp'));
% 
% fh_ARdpca = figure('position',[50,500,1000,600]); 
% hold on; 
% for icp=1:10
%     subplot(5,2,icp); hold on; grid on; box on
%     title(['PC of AR coeffs ' num2str(icp) ' , VE: ' num2str(ARdVarExp(icp))])
%     
%     plot(ARdPCs(:,icp)); xlim([0 arp])
% end
% set(fh_ARdpca,'color','w')

% AR_fh = figure('position',[50,500,1200,500]);
% imagesc(dARp')
% hold on; 
% title('AR(20) of deterended voxel time series.')
% colorbar; xlabel('voxels'); ylabel('AR coefficients')
% set(AR_fh,'color','w')

%%%% SAVE FIGURES %%%%%%%%%%%%%%%%%%%%%%%
% export_fig(fh_pca,'RFig/PC10.png')
% export_fig(fh_dpca,'RFig/dPC10.png')
% export_fig(AR_fh,'RFig/ARfig.png')
% export_fig(fh_ARdpca,'RFig/PCAAR20.png')
