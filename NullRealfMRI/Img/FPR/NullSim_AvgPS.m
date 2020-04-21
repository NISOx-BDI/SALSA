clear

Path2Mats='/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_ROCKLAND/'; 
FigDir = '/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/Img/FPR/RFigs';

load([Path2Mats '/ROCKLAND645SUBLIST.mat'])
SubList=ROCKLAND645SUBLIST(1:50); 
SubList(36) = []; 

for pwmethod = {'AR-W','AR-YW','ARMAHR','ACF'}
    psfh = figure('position',[50,500,1100,300]);

    AROlist = [1 2 5 10 20];
    MAOrd = 0;
    if strcmp(pwmethod{1},'ARMAHR'); MAOrd = 1; end; 
    if strcmp(pwmethod{1},'ACF'); AROlist = [5 10 15 30 60]; end; 
    
    gmean = []; 
    dt_cnt = 1; 
    for dtmethod = {'dct10','poly4','hpf100'}
        
        subplot(1,3,dt_cnt); 
        hold on; grid on; axis tight; box on; 
        title([pwmethod{1} ', ' dtmethod{1} ])
        pwRES = []; nRES = [];
        for AROcnt = 1:5
            disp([pwmethod{1} ', ' num2str(AROlist(AROcnt)) ', ' dtmethod{1} ])
            s_cnt = 1; 
            for s = SubList'
                %disp(s{1})
                A=[Path2Mats 'R.PW/' dtmethod{1} '/mats/' s{1} '_ses-DS2_EDboxcar_20_' pwmethod{1} '_AR' num2str(AROlist(AROcnt)) '_MA' num2str(MAOrd) '_FWHM0_' dtmethod{1} '.mat'];
                
                if ~exist(A,'file'); disp(['doesnot exists: ' A]); continue; end; 
                
                B=load(A,'SPEC');
                pwRES = [pwRES B.SPEC.Y_pwRES];
                nRES  = [nRES B.SPEC.Y_RES];
                gmean = [gmean geomean(B.SPEC.Y_pwRES(2:end))];

                s_cnt = s_cnt + 1; 
            end

            mpwRES = mean(pwRES,2);
            mnRES  = mean(nRES,2);
            plot(B.SPEC.X_RES(1:300),mpwRES(1:300),'LineWidth',1.3);
            ylim([0.01 0.04])
        end
        
        plot(B.SPEC.X_RES(1:300),mnRES(1:300),'LineWidth',1.3)
        plot(B.SPEC.X_RES(1:300),ones(1,300).*mean(gmean),'LineWidth',1.3,'color','k','linestyle','-.')
        
        legend([cellfun(@num2str,num2cell(AROlist),'UniformOutput',false),'Naive','geomean'])
        
        ylabel('Amplitude')
        xlabel('Frequency')
        
        
        dt_cnt = dt_cnt + 1; 
    end
    set(psfh,'color','w')
    export_fig(psfh,[FigDir '/Fig_AvgSpec_' pwmethod{1} '_FWHM0.png'])
end