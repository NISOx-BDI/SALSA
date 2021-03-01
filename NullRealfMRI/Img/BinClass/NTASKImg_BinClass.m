clear


% CohortID = 'ROCKLAND'; %ROCKLAND Cambridge Beijing
% T       = 900; 
% TR      = 0.645;

CohortID = 'NEO'; %ROCKLAND Cambridge Beijing
T       = 2300; 
TR      = 0.392;
SubNum  = 50; 

taskSNR = 'kernel'; %num2str(taskSNR.*100)
icaflag = 0;
gsrflag = 0;
FWHMl   = 5;
EDType  = 'ERF'; 
BCl     = 10;


Path2Mats=['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_' CohortID '/']; 
FigDir = '/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/Img/BinClass/RFigs11';
SubIDPath = [Path2Mats '/' CohortID '_subid.mat'];
load(SubIDPath)
SubList=participants(1:51); 

SesIDPath = [Path2Mats '/' CohortID '_sesid.mat'];
load(SesIDPath)
SesList=ses(1:SubNum);

Col = get(groot,'defaultAxesColorOrder');

% SubIDPath=['/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/Img/DL/' CohortID];
% load([SubIDPath '/participants.mat'])
% SubList=participants(1:51); 

pwmethodlist = {'ACF','gFAST','gACFadjxACFadj2tT1S0P5','AR-W','3dREMLfit'}; %{'ACF','ACFadj','gFAST','ACFadjx2','AR-W'};
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

elseif strcmp(CohortID,'NEO')
    dtla=[15,7,100,90];
    ACFOList = [5 10 15 fix(sqrt(T)) 2.*fix(sqrt(T)) -2];     
    
elseif strcmp(CohortID,'Beijing')
    dtla=[8,4,100,90]; 
    ACFOList = [5 10 fix(sqrt(T)) fix(2.*sqrt(T)) -2]; 
end

if gsrflag || icaflag || ~strcmp(CohortID,'ROCKLAND')
    GSRARMAflag = ['_GSR' num2str(gsrflag) '_AROMA' num2str(icaflag)];
else
    GSRARMAflag = '';
end

pwm_cnt = 1;
for pwmethod = pwmethodlist % {'AR-W','ACF','ACFadj'};%,'ARMAHR'}
    AROlist = [1 2 5 10 20];
    MAOrd = 0;

    if strcmp(pwmethod{1},'ARMAHR')
        MAOrd = 1; 
        AROlist = ARMAHRList;
    end; 
    
    if contains(pwmethod{1},'ACF')
        AROlist = ACFOList; 
    elseif strcmp(pwmethod{1},'gFAST')
        AROlist = 1;         
    elseif strcmp(pwmethod{1},'3dREMLfit')
        AROlist = 1;
        MAOrd = 1;
    end; 
    
    dt_cnt = 1; 
    for dtmethod = dtmethodlist
        for AROcnt = 1:numel(AROlist)
            disp([pwmethod{1} ', ' num2str(AROlist(AROcnt)) ', ' dtmethod{1} ])
            s_cnt = 1; 
            for s = 1:SubNum %s = SubList'
                %disp(s{1})
                %A=[Path2Mats 'R.PW/T' num2str(T) '/' dtmethod{1} '/sub-' s{1} '_ses-DS2_EDboxcar_20_' pwmethod{1} '_AR' num2str(AROlist(AROcnt)) '_MA' num2str(MAOrd) '_FWHM0_' dtmethod{1} '_tVALUE_PW' GSRARMAflag '.nii.gz'];
                
                A=[Path2Mats 'R.PW/T' num2str(T) '/' CohortID '_NTASK_' num2str(TR*1000) '_' num2str(T) '_' pwmethod{1} '_AR-' num2str(AROlist(AROcnt)) '_MA-' num2str(MAOrd) '_FWHM' num2str(FWHMl) '_' dtmethod{1} '_gsr' num2str(gsrflag) '_aroma' num2str(icaflag) '/' SubList{s} '_' SesList{s}];
                %EDERF_0_ACF_AR20_MA0_FWHM5_poly4_wtv_ICACLEAN1_GSR0.nii
                A=[A '/NTASK_SNR_' taskSNR '_ED'  EDType '_' num2str(BCl) '_' pwmethod{1} '_AR' num2str(AROlist(AROcnt)) '_MA' num2str(MAOrd) '_FWHM' num2str(FWHMl) '_' dtmethod{1} num2str(dtla(dt_cnt)) '_ICACLEAN' num2str(icaflag) '_GSR' num2str(gsrflag) '.mat'];
                              
                if ~exist(A,'file'); disp(['doesnot exists: ' A]); continue; end; elig{dt_cnt}{pwm_cnt}{AROcnt}{s_cnt} = s;
                load(A,'NTASK');
                %BIN{dt_cnt}{pwm_cnt}{AROcnt}(:,s_cnt)  = NTASK.BIN;
                wBIN{dt_cnt}{pwm_cnt}{AROcnt}(:,s_cnt) = NTASK.wBIN;
                s_cnt = s_cnt + 1; 
            end            
        end
        dt_cnt = dt_cnt + 1; 
    end
    pwm_cnt = pwm_cnt + 1; 
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


% dtmethodlist = {'dct','poly','hpf','hpfk'};
senspc_list = {'Sensitivity','Specificity','Accuracy'}; % 1: sensitivity, 2: specificity
for senspc_cnt = 1:numel(senspc_list)
        
    psfh = figure('position',[50,500,230,230]);
    hold on; grid on; box on; 
    title(['FWHM=' num2str(FWHMl) 'mm'], 'Interpreter','latex')
    
    fMRIStat = wBIN{2}{find(strcmpi(pwmethodlist,'AR-W'))}{1};
    FAST     = wBIN{1}{find(strcmpi(pwmethodlist,'gFAST'))}{1};
    FILM     = wBIN{3}{find(strcmpi(pwmethodlist,'ACF'))}{4};
    SALSA    = wBIN{2}{find(strcmpi(pwmethodlist,'gACFadjxACFadj2tT1S0P5'))}{6};
    REMLfit  = wBIN{2}{find(strcmpi(pwmethodlist,'3dREMLfit'))}{1};
    
    tmpplot = [mean(fMRIStat(senspc_cnt,:),2),mean(FAST(senspc_cnt,:),2),mean(FILM(senspc_cnt,:),2),mean(REMLfit(senspc_cnt,:),2),mean(SALSA(senspc_cnt,:),2)];
    bh0 = bar(tmpplot,'FaceColor',[.5 .5 .5]);
   
    ylabel(senspc_list{senspc_cnt},'FontSize',12,'Interpreter','latex')
    xlabel('Method','FontSize',12,'Interpreter','latex')
    psfh.Children.XTick=1:5;
    psfh.Children.XTickLabel={'fMRIStat','FAST','FILM','3dREMLfit','SALSA'};
    psfh.Children.XTickLabelRotation = 45;     


    bh0.FaceColor = Col(senspc_cnt,:);
    
    ctr = []; ydt = [];
    for k1 = 1:numel(pwmethodlist)
        ctr_tmp = bh0.XData;
        ydt_tmp = bh0.YData;

        ctr(k1,:) = ctr_tmp;
        ydt(k1,:) = ydt_tmp;

        for m_cnt = 1:numel(pwmethodlist)
            txth = text(ctr_tmp(m_cnt),ydt_tmp(m_cnt)+0.05,num2str(round(ydt_tmp(m_cnt),3)),'fontsize',10);
            set(txth,'rotation',90)
        end

    end    
    ylim([0 0.45]);

    set(psfh,'color','w')

    export_fig(psfh,[FigDir '/Fig_fprplots_ED' EDType '_' CohortID '_T' num2str(T) '_TR' num2str(TR.*1000) '_' senspc_list{senspc_cnt} '_FWHM' num2str(FWHMl) '_GSR' num2str(gsrflag) '_aroma' num2str(icaflag) '.png'])


end





% senspc_list = {'Sensitivity','Specificity','Accuracy'}; % 1: sensitivity, 2: specificity
% for senspc_cnt = 1:numel(senspc_list)
%         
%     psfh = figure('position',[50,500,1350,600]);
%     pwm_cnt = 1;
%     for pwmethod = pwmethodlist %{'AR-W','ACF','ACFadj'}%,'ARMAHR'}
% 
%         AROlist = [1 2 5 10 20];
%         MAOrd = 0;
%         
%         if strcmp(pwmethod{1},'ARMAHR')
%             MAOrd   = 1; 
%             AROlist = ARMAHRList;
%         end
%         
%         if contains(pwmethod{1},'ACF')
%             AROlist = ACFOList;  
%         elseif strcmp(pwmethod{1},'gFAST')
%             AROlist = 1;  
%         elseif strcmp(pwmethod{1},'3dREMLfit')
%             AROlist = 1;  
%             
%         end; 
% 
%         dt_cnt = 1; 
%         for dtmethod = dtmethodlist
%             for AROcnt = 1:numel(AROlist)
%                 disp([pwmethod{1} ', ' num2str(AROlist(AROcnt)) ', ' dtmethod{1} ])
%                 
%                 swBIN                         = wBIN{dt_cnt}{pwm_cnt}{AROcnt}; 
%                 saFPR(AROcnt)                 = mean(swBIN(senspc_cnt,:),2);
%             end
%             
%             sadFPR(:,dt_cnt)  = saFPR;
%             
%             clear saFPR
%             
%             dt_cnt = dt_cnt + 1; 
%         end
% 
%         spfh=subplot(2,ceil(numel(pwmethodlist)/2),pwm_cnt); 
%         hold on; grid on; box on; 
%         title(['  ' pwmethod{1} '$\%$, GSR' num2str(gsrflag) ', AROMA' num2str(icaflag) ', FWHMl' num2str(FWHMl)  ],'Interpreter','latex')
%         bar((sadFPR)) 
% 
%         if strcmp(pwmethod{1},'gFAST')
%             spfh.XTick=1:numel(dtmethodlist);
%             spfh.XTickLabel=dtmethodlist;
%             sph.XTickLabelRotation=45;
%             %line([0.5 numel(dtmethodlist)+.5],[alphlev alphlev],'LineWidth',1.3,'color','r')
%             xlabel('Detrending Method','Interpreter','latex')
%         else
%             spfh.XTick=1:numel(AROlist);
%             spfh.XTickLabel=cellfun(@num2str,num2cell(AROlist),'UniformOutput',false);
%             llgd = legend(dtmethodlist,'Location','southeast');
%             %line([0.5 numel(AROlist)+.5],[alphlev alphlev],'LineWidth',1.3,'color','r')
%             xlim([0.5 numel(AROlist)+.5])
%             xlabel('Model Order','Interpreter','latex')
%         end
%         
%         ylim([0 1])
%         ylabel(senspc_list{senspc_cnt},'Interpreter','latex')
% 
%         clear sadFPR
%         %clear sadFPR
%         pwm_cnt = pwm_cnt + 1; 
%     end
%     set(psfh,'color','w')
%     %export_fig(psfh,[FigDir '/Fig_fprplots_ED' EDType '_' CohortID '_T' num2str(T) '_TR' num2str(TR.*1000) '_' senspc_list{senspc_cnt} '_FWHM' num2str(FWHMl) '_GSR' num2str(gsrflag) '_aroma' num2str(icaflag) '.png'])
%     
% end