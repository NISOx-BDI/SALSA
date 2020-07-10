clear


clear

CohortID = 'Beijing'; %ROCKLAND Cambridge Beijing
T       = 225; 
TR      = 2;
icaflag = 0;
gsrflag = 0;
FWHMl   = 5;
EDType  = 'ER'; BCl = 0;


Path2Mats=['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_' CohortID '/']; 
FigDir = '/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/Img/FPR/RFigs3';
SubIDPath = [Path2Mats '/' CohortID '_subid.mat'];
load(SubIDPath)
SubList=participants(1:51); 

Col = get(groot,'defaultAxesColorOrder');

% SubIDPath=['/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/Img/DL/' CohortID];
% load([SubIDPath '/participants.mat'])
% SubList=participants(1:51); 

pwmethodlist = {'ACF','AR-W','gReMLxARw2t2j','gFAST','ACFadj','gFASTxACFadj','gFASTxACFadj2t','gFASTxACFadj2t2j'}; %{'ACF','ACFadj','gFAST','ACFadjx2','AR-W'};
dtmethodlist = {'dct','poly','hpf','hpfk'};

if strcmp(CohortID,'ROCKLAND')
    if T==900
        dtla=[10,4,100,90]; 
        ACFOList = [10 15 fix(sqrt(T)) fix(2.*sqrt(T)) -2]; 
    elseif T==404
        dtla=[9,4,100,90];  
        ACFOList = [10 15 fix(sqrt(T)) fix(2.*sqrt(T)) -2]; 
    end
elseif strcmp(CohortID,'Cambridge')
    dtla=[6,3,100];
    ACFOList = [5 10 20]; 
elseif strcmp(CohortID,'Beijing')
    dtla=[8,4,100,90]; 
    ACFOList = [5 10 fix(sqrt(T)) fix(2.*sqrt(T)) -2]; 
end

if gsrflag || icaflag || ~strcmp(CohortID,'ROCKLAND')
    GSRARMAflag = ['_GSR' num2str(gsrflag) '_AROMA' num2str(icaflag)];
else
    GSRARMAflag = '';
end

rankXlist   = [31,30,36,26,26,30];


get_tcrit = @(alph,df) tinv(1-alph,df);

pwm_cnt = 1;
for pwmethod = pwmethodlist % {'AR-W','ACF','ACFadj'};%,'ARMAHR'}
    AROlist = [1 2 5 10 20];
    MAOrd = 0;

    if strcmp(pwmethod{1},'ARMAHR')
        MAOrd = 1; 
        AROlist = ARMAHRList;
    end; 
    
    if strcmp(pwmethod{1},'ACF') || strcmp(pwmethod{1},'ACFadj') || strcmp(pwmethod{1},'ACFadjx2') || strcmp(pwmethod{1},'gFASTxACFadj') || strcmp(pwmethod{1},'gFASTxACFadj2t') || strcmp(pwmethod{1},'gFASTxACFadj2t2j')
        AROlist = ACFOList; 
    elseif strcmp(pwmethod{1},'gFAST')
        AROlist = 1;         
    end; 
    
    dt_cnt = 1; 
    for dtmethod = dtmethodlist
        for AROcnt = 1:numel(AROlist)
            disp([pwmethod{1} ', ' num2str(AROlist(AROcnt)) ', ' dtmethod{1} ])
            s_cnt = 1; 
            for s = SubList'
                %disp(s{1})
                %A=[Path2Mats 'R.PW/T' num2str(T) '/' dtmethod{1} '/sub-' s{1} '_ses-DS2_EDboxcar_20_' pwmethod{1} '_AR' num2str(AROlist(AROcnt)) '_MA' num2str(MAOrd) '_FWHM0_' dtmethod{1} '_tVALUE_PW' GSRARMAflag '.nii.gz'];
                
                A=[Path2Mats 'R.PW/T' num2str(T) '/' CohortID '_' num2str(TR*1000) '_' num2str(T) '_' pwmethod{1} '_AR-' num2str(AROlist(AROcnt)) '_MA-0_FWHM' num2str(FWHMl) '_' dtmethod{1} '_gsr' num2str(gsrflag) '_aroma' num2str(icaflag) '/' s{1} '_DS2'];
                %EDERF_0_ACF_AR20_MA0_FWHM5_poly4_wtv_ICACLEAN1_GSR0.nii
                A=[A '/ED'  EDType '_' num2str(BCl) '_' pwmethod{1} '_AR' num2str(AROlist(AROcnt)) '_MA0_FWHM' num2str(FWHMl) '_' dtmethod{1} num2str(dtla(dt_cnt)) '_wtv_ICACLEAN' num2str(icaflag) '_GSR' num2str(gsrflag) '.nii.gz'];
                              
                if ~exist(A,'file'); disp(['doesnot exists: ' A]); continue; end; elig{dt_cnt}{pwm_cnt}{AROcnt}{s_cnt} = s;
                tstat_tmp = CleanNIFTI_spm(A,'verbose',0);
                alltstat{dt_cnt}{pwm_cnt}{AROcnt}{s_cnt} = tstat_tmp;
                s_cnt = s_cnt + 1; 
            end            
        end
        dt_cnt = dt_cnt + 1; 
    end
    pwm_cnt = pwm_cnt + 1; 
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

alphlevlist = [0.05 0.01 0.001];
for alphlev_cnt = 1:numel(alphlevlist)
    
    alphlev = alphlevlist(alphlev_cnt);
    
    psfh = figure('position',[50,500,1350,600]);
    pwm_cnt = 1;
    for pwmethod = pwmethodlist %{'AR-W','ACF','ACFadj'}%,'ARMAHR'}

        AROlist = [1 2 5 10 20];
        MAOrd = 0;
        
        if strcmp(pwmethod{1},'ARMAHR')
            MAOrd   = 1; 
            AROlist = ARMAHRList;
        end
        
        if strcmp(pwmethod{1},'ACF') || strcmp(pwmethod{1},'ACFadj') || strcmp(pwmethod{1},'ACFadjx2') || strcmp(pwmethod{1},'gFASTxACFadj') || strcmp(pwmethod{1},'gFASTxACFadj2t') || strcmp(pwmethod{1},'gFASTxACFadj2t2j')
            AROlist = ACFOList;  
        elseif strcmp(pwmethod{1},'gFAST')
            AROlist = 1;  
        end; 

        dt_cnt = 1; 
        for dtmethod = dtmethodlist

    %         subplot(1,3,dt_cnt); 
    %         hold on; grid on; axis tight; box on; 
    %         title([pwmethod{1} ', ' dtmethod{1} ])
                df     = T-rankXlist(dt_cnt);
                t_crit = get_tcrit(alphlev,df);
            for AROcnt = 1:numel(AROlist)
                disp([pwmethod{1} ', ' num2str(AROlist(AROcnt)) ', ' dtmethod{1} ])
                
                num_sigvox = 0; allvox= 0 ; 
                for s_cnt = 1:numel(alltstat{dt_cnt}{pwm_cnt}{AROcnt})                    
                    tstat          = alltstat{dt_cnt}{pwm_cnt}{AROcnt}{s_cnt};
                    num_sigvox     = num_sigvox+sum(tstat>t_crit);
                    allvox         = allvox + numel(tstat); 
                    %sFPR(s_cnt)    = num_sigvox./numel(tstat);
                end
                
                sFPR                       = num_sigvox./allvox;
                FPR{pwm_cnt,AROcnt,dt_cnt} = sFPR;
                saFPR(AROcnt)              = sFPR;
            end
            
            sadFPR(:,dt_cnt)  = saFPR;
            
            clear saFPR
            
            dt_cnt = dt_cnt + 1; 
        end

        spfh=subplot(2,ceil(numel(pwmethodlist)/2),pwm_cnt); 
        hold on; grid on; box on; 
        title(['  ' pwmethod{1} ', $\alpha$=' num2str(alphlev*100) '$\%$, GSR' num2str(gsrflag) ', AROMA' num2str(icaflag) ', FWHMl' num2str(FWHMl)  ],'Interpreter','latex')
        bar(sadFPR) 

        if strcmp(pwmethod{1},'gFAST')
            spfh.XTick=1:numel(dtmethodlist);
            spfh.XTickLabel=dtmethodlist;
            sph.XTickLabelRotation=45;
            line([0.5 numel(dtmethodlist)+.5],[alphlev alphlev],'LineWidth',1.3,'color','r')
            xlabel('Detrending Method','Interpreter','latex')
        else
            spfh.XTick=1:numel(AROlist);
            spfh.XTickLabel=cellfun(@num2str,num2cell(AROlist),'UniformOutput',false);
            llgd = legend(dtmethodlist,'Location','northeast');
            line([0.5 numel(AROlist)+.5],[alphlev alphlev],'LineWidth',1.3,'color','r')
            xlim([0.5 numel(AROlist)+.5])
            xlabel('Model Order','Interpreter','latex')
        end
        
        ylim([0 alphlev.*2.*alphlev_cnt])
        ylabel('FPR','Interpreter','latex')

        clear sadFPR
        %clear sadFPR
        pwm_cnt = pwm_cnt + 1; 
    end
    set(psfh,'color','w')
    export_fig(psfh,[FigDir '/Fig_fprplots_ED' EDType '_' CohortID '_T' num2str(T) '_TR' num2str(TR.*1000) '_' num2str(alphlev.*1000) '_FWHM' num2str(FWHMl) '_GSR' num2str(gsrflag) '_aroma' num2str(icaflag) '.png'])
    
end