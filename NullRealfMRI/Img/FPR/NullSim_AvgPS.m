clear

CohortID = 'Cambridge'; %ROCKLAND Cambridge Beijing
T = 119;
icaflag = 1;
gsrflag = 1;

XLim    = 0.17; 

Path2Mats=['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_' CohortID '/']; 
FigDir = '/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/Img/FPR/RFigs';
SubIDPath = [Path2Mats '/' CohortID '_subid.mat'];
load(SubIDPath)
SubList=participants(1:51); 

Col = get(groot,'defaultAxesColorOrder');

% SubIDPath=['/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/Img/DL/' CohortID];
% load([SubIDPath '/participants.mat'])
% SubList=participants(1:51); 

if strcmp(CohortID,'ROCKLAND')
    dtmethodlist = {'dct10','poly4','hpf100'};
elseif strcmp(CohortID,'Cambridge')
    dtmethodlist = {'dct6','poly3','hpf100'};
    ACFOList = [5 10 20]; 
elseif strcmp(CohortID,'Beijing')
    dtmethodlist = {'dct8','poly4','hpf100'};
    ACFOList = [5 10 15 30]; 
end

if gsrflag || icaflag || ~strcmp(CohortID,'ROCKLAND')
    GSRARMAflag = ['_GSR' num2str(gsrflag) '_AROMA' num2str(icaflag)];
else
    GSRARMAflag = '';
end

for pwmethod = {'ACFadj','AR-W','ACF','ARMAHR'}
    psfh = figure('position',[50,500,1100,300]);

    AROlist = [1 2 5 10 20];
    MAOrd = 0;
    if strcmp(pwmethod{1},'ARMAHR'); MAOrd = 1; end; 
    if strcmp(pwmethod{1},'ACF') || strcmp(pwmethod{1},'ACFadj')
        AROlist = ACFOList;
    end 
    
    gmean = []; 
    dt_cnt = 1; 
    for dtmethod = dtmethodlist %{'dct10','poly4','hpf100'}
        
        subplot(1,3,dt_cnt); 
        hold on; grid on; axis tight; box on; 
        title([ CohortID ', T' num2str(T) ', ' pwmethod{1} ', ' dtmethod{1} ', GSR' num2str(gsrflag) ', AROMA' num2str(icaflag) ])
        pwRES = []; nRES = [];
        for AROcnt = 1:numel(AROlist)
            disp([pwmethod{1} ', ' num2str(AROlist(AROcnt)) ', ' dtmethod{1} ])
            s_cnt = 1; 
            for s = SubList'
                %disp(s{1})
                A=[Path2Mats 'R.PW/T' num2str(T) '/' dtmethod{1} '/mats/sub-' s{1} '_ses-DS2_EDboxcar_20_' pwmethod{1} '_AR' num2str(AROlist(AROcnt)) '_MA' num2str(MAOrd) '_FWHM0_' dtmethod{1} GSRARMAflag '.mat'];
                
                if ~exist(A,'file'); disp(['doesnot exists: ' A]); continue; end; 
                
                B=load(A,'SPEC');
                
                if numel(B.SPEC.Y_pwRES)<T/2+1; disp('too short.'); continue; end
                    
                pwRES = [pwRES B.SPEC.Y_pwRES(1:round(T/2))];
                nRES  = [nRES B.SPEC.Y_RES(1:round(T/2))];
                gmean = [gmean geomean(B.SPEC.Y_pwRES(2:end))];

                s_cnt = s_cnt + 1; 
            end

            mpwRES = mean(pwRES,2);
            mnRES  = mean(nRES,2);
            plot(B.SPEC.X_RES(1:round(T/2)),mpwRES,'LineWidth',1.3);
            ylim([0 1.5])
        end
        
        plot(B.SPEC.X_RES(1:round(T/2)),mnRES,'LineWidth',1.3,'color',Col(6,:))
        plot(B.SPEC.X_RES(1:round(T/2)),ones(length(B.SPEC.X_RES(1:round(T/2))),1),'color','k','linestyle','-.')
        %plot(B.SPEC.X_RES(1:400),ones.*mean(gmean),'LineWidth',1.3,'color','k','linestyle','-.')
        
        legend([cellfun(@num2str,num2cell(AROlist),'UniformOutput',false),'Naive','Expected SPD of GWN'],'Location','southeast')
        
        ylabel('Power')
        xlabel('Frequency')
        
        %xlim([0 0.1])
        xlim([0 XLim])
        
        dt_cnt = dt_cnt + 1; 
    end
    set(psfh,'color','w')
    export_fig(psfh,[FigDir '/Fig_AvgSpec_' CohortID '_'  num2str(T) '_' pwmethod{1} '_FWHM0' GSRARMAflag '_xlim' num2str(XLim*10)  '.pdf'])
end