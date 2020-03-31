clear

addpath('/Users/sorooshafyouni/Home/GitClone/DVARS/mis')
addpath('/Users/sorooshafyouni/Home/BCF/ARMA/Trend')
addpath('/Users/sorooshafyouni/Home/matlab/spm12')

SubID='A00029304';
SesID='DS2';

TR=0.645;

ImgDir='/Users/sorooshafyouni/Home/PREW/SimulatefMRI/R.mpp';
Path2Img=[ImgDir '/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-' num2str(TR*1000) '_bold.mpp/prefiltered_func_data_bet.nii.gz'];


%%% Read The Data %%%%%%%%%%%%%%%%%%%%%%%%
[Y,Stat]=CleanNIFTI_fsl(Path2Img,'demean','Norm',100);
T = Stat.CleanedDim(2);
V = Stat.CleanedDim(1);

Y = Y(4000:10000,:); 
V = size(Y,1);

%%% Generate a Design Matrix %%%%%%%%%%%%
bcduration  = 20;
BCWidth     = fix(bcduration/TR);
oneSTIM     = [zeros(1, BCWidth), ones(1, BCWidth)];
numberSTIM  = fix(T/numel(oneSTIM));
DD          = repmat(oneSTIM, [1, numberSTIM]);
hrf         = spm_hrf(TR); 
conv_dsgn   = conv(DD,hrf); 
X           = conv_dsgn(1:T)';
X           = X - mean(X); 

%%% DETREND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pp = 1+fix(TR.*T/150);
dY = multpolyfit(repmat(1:T,V,1),Y,T,pp);
dY = dY - mean(dY,2); 

%%% RESIDUALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
pinvX = pinv(X);
ResidFormingMat = eye(T)-X*pinvX;

residY = zeros(V,T);
for vi = 1:V    
    if ~mod(vi,10000); disp(['Residuals ::: on voxel ' num2str(vi)]); end;
    residY(vi,:) = ResidFormingMat*dY(vi,:)';
end

% Get ACF of residuals
ARord      = 20; 

disp(['Calculate the autocorrelation coefficients.'])
[~,~,dRESacov]     = AC_fft(residY,T);    

% Calculate the Autocorrelation Length
dRESacorr = dRESacov./sum(abs(residY).^2,2);
ACL = sum(dRESacorr.^2,2); 

YWflag = 0;
WrosleyFlag = 0;
ACFflag = 1;


%%% BIAS REDUCTION OF AUTOCORRELATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if WrosleyFlag
    M_biasred=zeros(ARord+1);
    for i=1:(ARord+1)
        Di=(diag(ones(1,T-i+1),i-1)+diag(ones(1,T-i+1),-i+1))/(1+(i==1));
        for j=1:(ARord+1)
           Dj=(diag(ones(1,T-j+1),j-1)+diag(ones(1,T-j+1),-j+1))/(1+(j==1));
           M_biasred(i,j)=trace(ResidFormingMat*Di*ResidFormingMat*Dj)/(1+(i>1));
        end
    end
    invM_biasred = inv(M_biasred);
else
    invM_biasred = eye(ARord+1);
end


YWARparam   = zeros(V,ARord);
dpwRES      = zeros(V,T); 
Bhatdpw     = zeros(V,1);

for vi = 1:V
    if ~mod(vi,1000); disp(['AR-YW ::: on voxel ' num2str(vi)]); end; 
    
    if YWflag 
        % Prewhiten the residuals using Yule-Walker estimates
        
        R   = toeplitz(dRESacorr(vi,1:ARord));
        r   = dRESacorr(vi,2:ARord+1);    
        YWARparam_tmp = R\r';
        YWARparam(vi,:) = YWARparam_tmp;
        ACMat     = full(spm_Q(YWARparam_tmp,T));
        
    elseif ACFflag 
        MTukey      = round(2*sqrt(T));
        %acfdRES     = dRESacov(vi,:)./sum(abs(residY(vi,:)).^2);
        acfdRES     = dRESacorr(vi,:); 
        acfdRES_tt  = [1 TukeyTaperMe(acfdRES(2:end),T-1,MTukey)];
        ACMat       = toeplitz(acfdRES_tt);  
        
    elseif WrosleyFlag
        
        dRESac_adj      = (invM_biasred*dRESacov)./sum(abs(residY(vi,:)).^2);
        [Ainvt posdef]  = chol(toeplitz(dRESac_adj));
        nl      = size(Ainvt,1);
        A       = inv(Ainvt');
        B       = ones(n-nl,1)*A(nl,:);
        Vmhalf  = spdiags(B,1:nl,T-nl,T);
    end
    
    invACMat        = inv(ACMat); % pinv is damn slow!
    sqrtmpinvACMat  = chol(invACMat); % also sqrtm is really slow!
    %sqrtmpinvACMat = sqrtm(invACMat);
    
    % Make the X & Y whitened and re-estimate the Betas
    Ystar_YW = sqrtmpinvACMat*dY(vi,:)';
    Xstar_YW = sqrtmpinvACMat*X;
    
    % Fit a model to the prewhitened system
    BhatStar_tmp = Xstar_YW\Ystar_YW; % \Beta = X^{+}Y -- sanitycheck: [b,c,d] = glmfit(Xdpw,YYYdpw);
    YhatStar     = Xstar_YW*BhatStar_tmp; %
    dpwRES_tmp   = Ystar_YW-YhatStar; % sanitycheck: d.resid(1:5)
    
    Beta_PW_SE_S(vi) = BhatStar_tmp;
    dpwRES(vi,:)     = dpwRES_tmp;
    
    % estimate the theoritical se of the system
    Beta_PW_SE_T(vi) = sqrt(((dpwRES_tmp'*dpwRES_tmp))./(X'*X)/(T-2));
        
    % Whithout prewhitening of error:
    Bhat_tmp = X\dY(vi,:)'; % \Beta = X^{+}Y -- sanitycheck: [b,c,d] = glmfit(Xdpw,YYYdpw);
    Yhat = X*Bhat_tmp; %
    res  = dY(vi,:)'-Yhat; % sanitycheck: d.resid(1:5)
    Bhat_RAW_S(vi)    = Bhat_tmp; 
    Beta_RAW_SE_T(vi) = sqrt(((res'*res))./(X'*X)/(T-2));
        
end

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
