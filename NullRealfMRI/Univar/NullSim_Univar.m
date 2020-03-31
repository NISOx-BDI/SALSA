clear

addpath('/Users/sorooshafyouni/Home/GitClone/DVARS/mis')
addpath('/Users/sorooshafyouni/Home/BCF/ARMA/Trend')
addpath('/Users/sorooshafyouni/Home/matlab/spm12')
addpath('/Users/sorooshafyouni/Home/BCF/ARMA/AR_YW')

pwdmethod = 'AR-YW'; %ACF AR-YW AR-W

% Order of the AR which will be fitted to the data. 
Mord      = 5; 


rngid = 20200330;
rng(rngid)

SubID='A00029304';
SesID='DS2';
TR=0.645;

ImgDir='/Users/sorooshafyouni/Home/PREW/SimulatefMRI/R.mpp';
Path2Img=[ImgDir '/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-' num2str(TR*1000) '_bold.mpp/prefiltered_func_data_bet.nii.gz'];


%%% Read The Data %%%%%%%%%%%%%%%%%%%%%%%%
[Y,Stat]=CleanNIFTI_fsl(Path2Img,'demean');
T = Stat.CleanedDim(2);
Vorig = Stat.CleanedDim(1);

Y = Y(32543,:);  
Y = Y - mean(Y);
Y = Y./std(Y); 

Vorig = size(Y,1);

%%% Generate a Design Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bcduration  = 20;
BCWidth     = fix(bcduration/TR);
oneSTIM     = [zeros(1, BCWidth), ones(1, BCWidth)];
numberSTIM  = fix(T/numel(oneSTIM));
DD          = repmat(oneSTIM, [1, numberSTIM]);
hrf         = spm_hrf(TR); 
conv_dsgn   = conv(DD,hrf); 
X           = conv_dsgn(1:T)';
X           = X - mean(X); 

%%% DETREND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pp = 1+fix(TR.*T/150);
dYorig = multpolyfit(repmat(1:T,Vorig,1),Y,T,pp);
dYorig = dYorig - mean(dYorig,2); 

%%% SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Order of AR which will be used to simulated stuff
ARporder = 20;
ARParam  = AR_YW(dYorig,ARporder);
nRlz = 1000; 
model = arima('Constant',0,'AR',ARParam,'Variance',1);
dY = simulate(model,T,'NumPaths',nRlz)'; 
dY = dY - mean(dY,2); 

%%% RESIDUALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pinvX = pinv(X);
ResidFormingMat = eye(T)-X*pinvX;

residY = zeros(nRlz,T);
for vi = 1:nRlz    
    if ~mod(vi,10000); disp(['Residuals ::: on voxel ' num2str(vi)]); end;
    residY(vi,:) = ResidFormingMat*dY(vi,:)';
end

%%% ACFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Calculate the autocorrelation coefficients.'])
[~,~,dRESacov]     = AC_fft(residY,T);    

% Calculate the Autocorrelation Length
dRESacorr = dRESacov./sum(abs(residY).^2,2);
ACL = sum(dRESacorr.^2,2); 

%%% BIAS REDUCTION OF AUTOREGRESSIVE MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YWflag      = 0;
WrosleyFlag = 0;
ACFflag     = 0;
if strcmpi(pwdmethod,'AR-W')
    WrosleyFlag = 1; 
    warning('off','MATLAB:toeplitz:DiagonalConflict')
    M_biasred = zeros(Mord+1);
    for i=1:(Mord+1)
        Di = (diag(ones(1,T-i+1),i-1)+diag(ones(1,T-i+1),-i+1))/(1+(i==1));
        for j=1:(Mord+1)
           Dj = (diag(ones(1,T-j+1),j-1)+diag(ones(1,T-j+1),-j+1))/(1+(j==1));
           M_biasred(i,j) = trace( ResidFormingMat*Di*ResidFormingMat*Dj )/(1+(i>1));
        end
    end
    invM_biasred = inv(M_biasred);
elseif strcmpi(pwdmethod,'AR-YW')
    YWflag  = 1; 
    invM_biasred = eye(Mord+1);
elseif strcmpi(pwdmethod,'ACF') 
    ACFflag = 1;
end

for vi = 1:nRlz
    
    if YWflag 
        if ~mod(vi,500); disp(['AR-YW ::: on voxel ' num2str(vi)]); end; 
        % Prewhiten the residuals using Yule-Walker estimates
        YWARparam_tmp   = AR_YW(dYorig,Mord);
        YWARparam(vi,:) = YWARparam_tmp;
        ACMat           = full(spm_Q(YWARparam_tmp,T));
        
        invACMat        = inv(ACMat); % pinv is damn slow!
        %sqrtmpinvACMat  = chol(invACMat); % also sqrtm is really slow!
        sqrtmVhalf      = chol(invACMat);
        
    elseif ACFflag 
        if ~mod(vi,500); disp(['ACF-TUKEY ::: on voxel ' num2str(vi)]); end; 
        acfdRES     = dRESacorr(vi,:); 
        %Mord        = 2*round(sqrt(T));
        acfdRES_tt  = [1 TukeyTaperMe(acfdRES(2:end),T-1,Mord)];
        ACMat       = toeplitz(acfdRES_tt);  
        
        invACMat       = inv(ACMat); % pinv is damn slow!
        %sqrtmpinvACMat  = chol(invACMat); % also sqrtm is really slow!
        sqrtmVhalf = chol(invACMat);
        
    elseif WrosleyFlag
        if ~mod(vi,500); disp(['AR-Adj ::: on voxel ' num2str(vi)]); end; 
        dRESac_adj      = (invM_biasred*dRESacov(vi,1:Mord+1)');
        dRESac_adj      = dRESac_adj./dRESac_adj(1); % make a auto*correlation*
        
        [Ainvt,posdef] = chol(toeplitz(dRESac_adj)); 
        p1      = size(Ainvt,1); 
        A       = inv(Ainvt'); 
        sqrtmVhalf  = toeplitz([A(p1,p1:-1:1) zeros(1,T-p1)],zeros(1,T)); 
        sqrtmVhalf(1:p1,1:p1) = A;

    end
    
    % [inverse] covariance Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    
    % Make the X & Y whitened %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ystar_YW = sqrtmVhalf*dY(vi,:)';
    Xstar_YW = sqrtmVhalf*X;
    
    % Fit a model to the prewhitened system  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BhatStar_tmp = Xstar_YW\Ystar_YW; % \Beta = X^{+}Y 
    YhatStar     = Xstar_YW*BhatStar_tmp; %
    dpwRES_tmp   = Ystar_YW-YhatStar; % sanitycheck: d.resid(1:5)
    Beta_PW_SE_S(vi) = BhatStar_tmp;
    Beta_PW_SE_T(vi) = sqrt(((dpwRES_tmp'*dpwRES_tmp))./(Xstar_YW'*Xstar_YW)/(T-rank(X)));
    
            
    % Whithout prewhitening of error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Bhat_tmp    = X\dY(vi,:)'; % \Beta = X^{+}Y 
    Yhat        = X*Bhat_tmp; %
    resNaive    = dY(vi,:)'-Yhat; % sanitycheck: d.resid(1:5)
    Bhat_Naive_S(vi)    = Bhat_tmp; % SE
    Beta_Naive_SE_T(vi) = sqrt(((resNaive'*resNaive))./(X'*X)/(T-rank(X))); % SE
    
    [~,CPSstat_PW(vi),CPZ_PW(vi)] = CPSUnivar(dpwRES_tmp,X);
    [~,CPSstat_Naive(vi),CPZ_Naive(vi)] = CPSUnivar(resNaive,X);
        
end

Beta_PW_SE_S_STD    = std(Beta_PW_SE_S); 
Beta_PW_SE_T_MEAN   = mean(Beta_PW_SE_T);

Bhat_NAIVE_S_STD      = std(Bhat_Naive_S); 
Beta_NAIVE_SE_T_MEAN  = mean(Beta_Naive_SE_T); 


save(['R/NullSim_Univar_Biases_' pwdmethod '_AR' num2str(Mord) '.mat'],...
    'Beta_PW_SE_S','rngid','Beta_PW_SE_T','Bhat_Naive_S','Beta_Naive_SE_T',...
    'CPSstat_Naive','CPZ_Naive','CPSstat_PW','CPZ_PW',...
    'dY','Y','T','TR')

BetaBias = @(Truth_MC,Estimated_TH) ((Estimated_TH-Truth_MC)./Truth_MC).*100;
fh = figure; 
hold on; grid on; box on; 
title('[THEORETICAL - (TRUE) SIMULATION]/ (TRUE) SIMULATION x 100')
bh = bar([BetaBias(Bhat_NAIVE_S_STD,Beta_NAIVE_SE_T_MEAN),BetaBias(Beta_PW_SE_S_STD,Beta_PW_SE_T_MEAN)]);
fh.Children.XTick=1:2;
fh.Children.XTickLabel = {'Naive','Prewhitened'};
ylabel('Bias in SE of $\hat\beta$','Interpreter','latex')


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
