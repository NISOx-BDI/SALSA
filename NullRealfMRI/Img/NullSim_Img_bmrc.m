%clear

warning('off','all')

% ---------------- TEST ----------------------------
%pwdmethod  = 'ACF'; %ACF AR-YW AR-W ARMAHR
%Mord       = 30; 
%MPparamNum = 0;
%TempTreMethod = 'dct'
%SubID     = 'A00029304';
%SesID     = 'DS2';
%lFWHM     = 0;
%TR        = 0.645;
%COHORTDIR = '/well/nichols/users/scf915/ROCKLAND';
%Path2ImgResults = [COHORTDIR '/R.PW/TEST' pwdmethod '_AR-' num2str(Mord) '_MA-' num2str(MPparamNum)  ]

% What is flowing in from the cluster:
disp('From the cluster ======================')
disp(['SubID: ' SubID])
disp(['SesID: ' SesID])
disp(['TR: ' num2str(TR)])
disp(['ARmethod: ' pwdmethod])
disp(['AR order:' num2str(Mord)])
disp(['MA order: ' num2str(MPparamNum)])
disp(['lFWHM: ' num2str(lFWHM)])
disp(['Detrending: ' TempTreMethod])
disp(['COHORT directory:' COHORTDIR])
disp(['Parth 2 Results: ' Path2ImgResults ])
disp('=======================================')

PATH2AUX='~/bin/FILM2';
addpath([PATH2AUX '/utils/Trend'])
addpath([PATH2AUX '/utils/AR_YW'])
addpath([PATH2AUX '/utils/ARMA_HR'])
addpath([PATH2AUX '/mis'])
addpath (fullfile ('/users/nichols/scf915', 'spm12-r7771'));
%addpath('/well/nichols/users/scf915/externals/spm12')

disp('=====SET UP PATHS =============================')
%Raw Images (MMP feat output)
Path2ImgRaw = [COHORTDIR '/R_mpp'];
Path2ImgDir = [Path2ImgRaw '/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-' num2str(TR*1000) '_bold_mpp'];

if ~lFWHM
    Path2Img    = [Path2ImgDir '/prefiltered_func_data_bet.nii'];
else
    Path2Img    = [Path2ImgDir '/prefiltered_func_data_bet_FWHM' num2str(lFWHM) '.nii'];
end
  
Path2MC     = [Path2ImgDir '/prefiltered_func_data_mcf.par'];

disp(['Image: ' Path2Img])
disp(['Motion params: ' Path2MC])

% Directory 2 save the results
Path2ImgResults=[Path2ImgResults '/' SubID '_' SesID];
if ~exist(Path2ImgResults, 'dir')
	mkdir(Path2ImgResults)
	disp(['The directory: ' Path2ImgResults ' did not exists. I made one. '])
end

disp(['Output stuff: ' Path2ImgResults])

SaveImagesFlag      = 1; 
DoDetrendingPrior   = 0; 
MParamNum           = 24; 

%%% Read The Data %%%%%%%%%%%%%%%%%%%%%%%%
disp('=====LOAD THE IMAGE ===========================')
[Y,InputImgStat]=CleanNIFTI_spm(Path2Img,'demean');
T = InputImgStat.CleanedDim(2);
TR = InputImgStat.voxelsize(4);
Vorig = InputImgStat.CleanedDim(1);
V = Vorig;

if size(Y,1)~=T; Y = Y'; end; %TxV

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DESIGN MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++ Construct a design matrix')
%%% Generate a Design Matrix --------------------------------
EDtype = 'boxcar'; 
BCl = 20;
EDX = GenerateED(BCl,T,TR); 
EDX = EDX - mean(EDX); 

X   = EDX;
disp(['design updated, ' num2str(size(X,2))])

% Motion parameters ----------------------------------------
MCp = [];
if MParamNum     == 6
    MCp = load(Path2MC);
elseif MParamNum == 12
    MCp = load(Path2MC);
    MCp = [MCp,MCp.^2]; % 12 parameter motion 
elseif MParamNum == 24
    o6MCp   = load(Path2MC);
    so6MCp  = o6MCp.^2;
    do6MCp  = [o6MCp(1,:);  diff(o6MCp)]; % is this stupid?
    ds6MCp  = [so6MCp(1,:); diff(so6MCp)];% is this stupid?
    MCp     = [o6MCp,so6MCp,do6MCp,ds6MCp];
end
disp(['Number of motion parameter: ' num2str(MParamNum)])

X = [X,MCp];
disp(['design updated, ' num2str(size(X,2))])

% Temporal trends ----------------------------------------
TempTrend = [];
if ~exist('NumTmpTrend','var'); NumTmpTrend=[]; end;
if any(strcmpi(TempTreMethod,{'dct','spline','poly'}))
    [TempTrend,NumTmpTrend]   = GenerateTemporalTrends(T,TR,TempTreMethod,NumTmpTrend); % DC + Poly trends + Spline trends 
    TempTrend   = TempTrend(:,2:end); % we add a column of one later.
    
elseif strcmpi(TempTreMethod,{'hpf'})
    if isempty(NumTmpTrend) || ~exist('NumTmpTrend','var'); NumTmpTrend=100; end; 
    hp_ff = hp_fsl(T,NumTmpTrend,TR);
    
    X     = hp_ff*X;    % high pass filter the design
    dY    = hp_ff*dY;  % high pass filter the data
end
disp(['Detrending: ' TempTreMethod ',param: ' num2str(NumTmpTrend)])
%
X           = [X,TempTrend];
disp(['design updated, ' num2str(size(X,2))])

% Centre the design  ----------------------------------
X           = X - mean(X); % demean everything 

%%% RESIDUALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++Get the residuals using I-XX+')
pinvX           = pinv(X); 
ResidFormingMat = eye(T)-X*pinvX; % residual forming matrix 
residY          = ResidFormingMat*dY;

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
    YWflag          = 1; 
    invM_biasred    = eye(Mord+1);
elseif strcmpi(pwdmethod,'ACF') % ACF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ACFflag         = 1;
elseif strcmpi(pwdmethod,'ARMAHR') % ARMAHR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ARMAHRflag      = 1; 
    MPparamNum      = 1;  % the MA order 
    ARParamARMA     = 50; % the higher fit in ARMA HR
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIT A MODEL TO THE ACd DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++Fit data to the Naive model.')
X0                                          = [ones(T,1),X];
[Bhat_Naive,~,resNaive,Stat_Naive_SE_tmp]   = myOLS(dY,X0);
SE_Naive                                    = Stat_Naive_SE_tmp.se;
tVALUE_Naive                                = Stat_Naive_SE_tmp.tval;
[~,CPSstat_Naive,CPZ_Naive]                 = CPSUnivar(resNaive,X0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AUTOCORR & AUTOCOV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% ACFs %%%%%%%%%%%%%%%
disp(['++++++++++++Calculate the autocorrelation coefficients.'])
[~,~,dRESacov]  = AC_fft(residY,T); % Autocovariance; VxT
dRESacov        = dRESacov'; %TxV
dRESacorr       = dRESacov./sum(abs(residY).^2); % Autocorrelation
ACL             = sum(dRESacorr.^2); % Autocorrelation Length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PREWHITEN THE RESIDULAS & ESTIMATE BIAS AND CPS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Preallocate memory
Bhat_PW    = zeros(V,1);
SE_PW      = zeros(V,1); 
tVALUE_PW  = zeros(V,1);
CPSstat_PW = zeros(V,1); 
CPZ_PW     = zeros(V,1);

disp('++++++++++++Starts the voxel-wise prewhitening')

for vi = 1:V
    if YWflag % Yule-Walker %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp(['AR-YW ::: on voxel ' num2str(vi)]); end;        
        %AR_YW -------------------------------------------
        ac_tmp          = dRESacorr(:,vi);    
        R_tmp           = toeplitz(ac_tmp(1:Mord));
        r_tmp           = ac_tmp(2:Mord+1);
        %YWARparam_tmp   = pinv(R_tmp)*r_tmp;
        YWARparam_tmp   = R_tmp\r_tmp;
        % ------------------------------------------------
        
        ACMat           = full(spm_Q(YWARparam_tmp,T));
        
        invACMat        = inv(ACMat); % pinv is damn slow!
        sqrtmVhalf      = chol(invACMat);
        
    elseif ACFflag % ACF - Tukey tapered %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp(['ACF-TUKEY ::: on voxel ' num2str(vi)]); end; 
        acfdRES     = dRESacorr(:,vi); 
        %Mord        = 2*round(sqrt(T));
        acfdRES_tt  = [1 TukeyTaperMe(acfdRES(2:end),T-1,Mord)];
        ACMat       = toeplitz(acfdRES_tt);  
        
        invACMat   = inv(ACMat); % pinv is damn slow!
        sqrtmVhalf = chol(invACMat);
        
    elseif WrosleyFlag % Worsely %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp(['AR-W ::: on voxel ' num2str(vi)]); end; 
        dRESac_adj      = (invM_biasred*dRESacov(1:Mord+1,vi));
        dRESac_adj      = dRESac_adj./dRESac_adj(1); % make a auto*correlation*
        
        [Ainvt,posdef]      = chol(toeplitz(dRESac_adj)); 
        p1                  = size(Ainvt,1); 
        A                   = inv(Ainvt'); 
        sqrtmVhalf          = toeplitz([A(p1,p1:-1:1) zeros(1,T-p1)],zeros(1,T)); 
        sqrtmVhalf(1:p1,1:p1) = A;
        
    elseif ARMAHRflag % ARMA HR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~mod(vi,5000); disp(['ARMA-HR ::: on voxel ' num2str(vi)]); end; 
  
        %AR_YW -------------------------------------------
        ac_tmp              = dRESacorr(:,vi);    
        R_tmp               = toeplitz(ac_tmp(1:ARParamARMA));
        r_tmp               = ac_tmp(2:ARParamARMA+1);
        %YWARparam_tmp       = pinv(R_tmp)*r_tmp %
        YWARparam_tmp       = R_tmp\r_tmp;
        % ------------------------------------------------
        
        [arParam,maParam]   = ARMA_HR_ACF(residY(:,vi),YWARparam_tmp',T,Mord,MPparamNum);
        ACMat               = ARMACovMat([arParam,maParam],T,Mord,MPparamNum);

        invACMat   = inv(ACMat); % pinv is damn slow!
        sqrtmVhalf = chol(invACMat);        
    end
        
    % Make the X & Y whitened %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ystar_YW = sqrtmVhalf*dY(:,vi);
    Xstar_YW = sqrtmVhalf*X;   
    % Fit a model to the prewhitened system  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    Xstar_YW                                      = [ones(T,1), Xstar_YW]; % add intercept
    [Bhat_PW_S_tmp,~,dpwRES_tmp,Stat_PW_SE_T_tmp] = myOLS(Ystar_YW,Xstar_YW);
    Bhat_PW(vi)    = Bhat_PW_S_tmp;         
    SE_PW(vi)      = Stat_PW_SE_T_tmp.se;
    tVALUE_PW(vi)  = Stat_PW_SE_T_tmp.tval;
   
    [~,CPSstat_PW(vi),CPZ_PW(vi)] = CPSUnivar(dpwRES_tmp,X0);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE THE RESULTS AS AN IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++Save the results.')

if SaveImagesFlag
    % 3D IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VariableList = {'Bhat_Naive','SE_Naive','tVALUE_Naive',...
        'Bhat_PW','SE_PW','tVALUE_PW',...
        'CPSstat_PW','CPZ_PW',...
        'CPZ_Naive','CPSstat_Naive',...
        'ACL'};
    OutputImgStat            = InputImgStat.spmV(1);
    OutputImgStat.Removables = InputImgStat.Removables;

    for vname = VariableList

        tmpvar                   = eval(vname{1});
        OutputImgStat.fname      = [Path2ImgResults '/ED' EDtype '_' num2str(BCl) '_' pwdmethod '_AR' num2str(Mord) '_MA' num2str(MPparamNum) '_FWHM' num2str(lFWHM) '_' TempTreMethod num2str(NumTmpTrend) '_' vname{1} '.nii'];

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
