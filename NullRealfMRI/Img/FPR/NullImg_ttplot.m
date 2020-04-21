clear

Path2Mats='/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_ROCKLAND/'; 
FigDir = '/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/Img/FPR/RFigs';

load([Path2Mats '/ROCKLAND645SUBLIST.mat'])
SubList=ROCKLAND645SUBLIST(1:50); 
SubList(36) = []; 
rankXlist   = [31,30,36,26];
T = 900; 

get_tcrit = @(alph,df) tinv(1-alph,df);


pwm_cnt = 1;
for pwmethod = {'AR-W','AR-YW','ARMAHR','ACF'}
    AROlist = [1 2 5 10 20];
    MAOrd = 0;
    if strcmp(pwmethod{1},'ARMAHR'); MAOrd = 1; end; 
    if strcmp(pwmethod{1},'ACF'); AROlist = [5 10 15 30 60]; end; 
    dt_cnt = 1; 
    for dtmethod = {'dct10','poly4','hpf100'}
        for AROcnt = 1:5
            disp([pwmethod{1} ', ' num2str(AROlist(AROcnt)) ', ' dtmethod{1} ])
            s_cnt = 1; 
            for s = SubList'
                %disp(s{1})
                A=[Path2Mats 'R.PW/' dtmethod{1} '/' s{1} '_ses-DS2_EDboxcar_20_' pwmethod{1} '_AR' num2str(AROlist(AROcnt)) '_MA' num2str(MAOrd) '_FWHM0_' dtmethod{1} '_tVALUE_PW.nii.gz'];
                if ~exist(A,'file'); disp(['doesnot exists: ' A]); continue; end; 
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

%alphlevlist = [0.05 0.01 0.001];

%alphlev = alphlevlist(alphlev_cnt);
cvt_05 = get_tcrit(.05/2,T-30);
cvt_01 = get_tcrit(.01/2,T-30);
cvt_001 = get_tcrit(.001/2,T-30);

pwm_cnt = 1;
for pwmethod = {'AR-W','AR-YW','ARMAHR','ACF'}
    psfh = figure('position',[50,500,1200,700]);

    AROlist = [1 2 5 10 20];
    MAOrd = 0;
    if strcmp(pwmethod{1},'ARMAHR'); MAOrd = 1; end; 
    if strcmp(pwmethod{1},'ACF'); AROlist = [5 10 15 30 60]; end; 

    dt_cnt = 1; 
    for dtmethod = {'dct10','poly4','hpf100'}

%         subplot(1,3,dt_cnt); 
%         hold on; grid on; axis tight; box on; 
%         title([pwmethod{1} ', ' dtmethod{1} ])
        for AROcnt = 1:5
            disp([pwmethod{1} ', ' num2str(AROlist(AROcnt)) ', ' dtmethod{1} ])

            for s_cnt = 1:numel(alltstat{dt_cnt}{pwm_cnt}{AROcnt})
%                     A=[Path2Mats 'R.PW/' dtmethod{1} '/' s{1} '_ses-DS2_EDboxcar_20_' pwmethod{1} '_AR' num2str(AROlist(AROcnt)) '_MA' num2str(MAOrd) '_FWHM0_' dtmethod{1} '_tVALUE_PW.nii.gz'];
%                     if ~exist(A,'file'); disp(['doesnot exists: ' A]); continue; end; 

                tstat          = alltstat{dt_cnt}{pwm_cnt}{AROcnt}{s_cnt};
                df             = T-rankXlist(dt_cnt);

                [ttx(s_cnt,:),tty(s_cnt,:)]  = tqqplot(tstat,df);

            end

            sattx(:,AROcnt) = mean(ttx);
            satty(:,AROcnt) = mean(tty);
        end

        spfhRow1=subplot(2,3,dt_cnt);
        hold on; grid on; box on; 
        title([pwmethod{1} ', ' dtmethod{1} ', Right Tail'],'Interpreter','latex')
        plot(sattx,satty,'LineWidth',1.3)
        plot(0:4,0:4,'LineWidth',1.3,'LineStyle','-.','color',[.5 .5 .5])
        ylim([0 4])
        xlim([0 4])
        
        line([cvt_05 cvt_05],[0 4],'LineWidth',1.1,'LineStyle','-.','color','r')
        line([cvt_01 cvt_01],[0 4],'LineWidth',1.1,'LineStyle','-.','color','r')
        line([cvt_001 cvt_001],[0 4],'LineWidth',1.1,'LineStyle','-.','color','r')
        
        xlabel('tinv','Interpreter','latex')
        ylabel('t-statistics','Interpreter','latex')

        legend(cellfun(@num2str,num2cell(AROlist),'UniformOutput',false),'Location','northwest');
        
        spfhRow2=subplot(2,3,dt_cnt+3);
        hold on; grid on; box on; 
        title([pwmethod{1} ', ' dtmethod{1} ', Left Tail'],'Interpreter','latex')
        plot(sattx,satty,'LineWidth',1.3)
        plot(-4:0,-4:0,'LineWidth',1.3,'LineStyle','-.','color',[.5 .5 .5])
        ylim([-4 0])
        xlim([-4 0])
        
        line(-[cvt_05 cvt_05],[-4 0],'LineWidth',1.1,'LineStyle','-.','color','r')
        line(-[cvt_01 cvt_01],[-4 0],'LineWidth',1.1,'LineStyle','-.','color','r')
        line(-[cvt_001 cvt_001],[-4 0],'LineWidth',1.1,'LineStyle','-.','color','r')
        
        xlabel('tinv','Interpreter','latex')
        ylabel('t-statistics','Interpreter','latex')
        
        dt_cnt = dt_cnt + 1; 
    end

    %clear sadFPR
    pwm_cnt = pwm_cnt + 1; 
    set(psfh,'color','w')
    export_fig(psfh,[FigDir '/Fig_ttplots_' pwmethod{1} '_FWHM0.png'])
end
