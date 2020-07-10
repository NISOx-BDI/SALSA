clear

Path2Mats='/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_ROCKLAND/'; 
FigDir = '/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/Img/FPR/RFigs';

SubIDPath='/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/Img/DL/ROCKLAND';
load([SubIDPath '/participants.mat'])
SubList=participants(1:50); 

rankXlist   = [31,30,36,26];
T = 404; 
dtmethodlist = {'dct9','poly4','hpf100'};

get_tcrit = @(alph,df) tinv(1-alph,df);


pwm_cnt = 1;
for pwmethod = {'AR-W','ACFadj','ACF','ARMAHR'}
    AROlist = [1 2 5 10 20];
    MAOrd = 0;
    if strcmp(pwmethod{1},'ARMAHR'); MAOrd = 1; end; 
    if strcmp(pwmethod{1},'ACF') || strcmp(pwmethod{1},'ACFadj')
        AROlist = [5 10 15 fix(sqrt(T)) fix(2.*sqrt(T))]; 
    end 
    
    dt_cnt = 1; 
    for dtmethod = dtmethodlist
        for AROcnt = 1:5
            disp([pwmethod{1} ', ' num2str(AROlist(AROcnt)) ', ' dtmethod{1} ])
            s_cnt = 1; 
            for s = SubList'
                %disp(s{1})
                A=[Path2Mats 'R.PW/T' num2str(T) '/' dtmethod{1} '/sub-' s{1} '_ses-DS2_EDboxcar_20_' pwmethod{1} '_AR' num2str(AROlist(AROcnt)) '_MA' num2str(MAOrd) '_FWHM0_' dtmethod{1} '_tVALUE_PW.nii.gz'];
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
cvt_05  = -log10(.05);
cvt_01  = -log10(.01);
cvt_001 = -log10(.001);

pwm_cnt = 1;
for pwmethod = {'AR-W','ACFadj','ACF','ARMAHR'}
    psfh = figure('position',[50,500,1200,700]);

    AROlist = [1 2 5 10 20];
    MAOrd = 0;
    if strcmp(pwmethod{1},'ARMAHR'); MAOrd = 1; end; 
    if strcmp(pwmethod{1},'ACF') || strcmp(pwmethod{1},'ACFadj')
        AROlist = [5 10 15 fix(sqrt(T)) fix(2.*sqrt(T))]; 
    end 
    dt_cnt = 1; 
    for dtmethod = dtmethodlist

%         subplot(1,3,dt_cnt); 
%         hold on; grid on; axis tight; box on; 
%         title([pwmethod{1} ', ' dtmethod{1} ])
        for AROcnt = 1:5
            disp([ 'T' num2str(T) ', ' pwmethod{1} ', ' num2str(AROlist(AROcnt)) ', ' dtmethod{1} ])

            for s_cnt = 1:numel(alltstat{dt_cnt}{pwm_cnt}{AROcnt})
%                     A=[Path2Mats 'R.PW/' dtmethod{1} '/' s{1} '_ses-DS2_EDboxcar_20_' pwmethod{1} '_AR' num2str(AROlist(AROcnt)) '_MA' num2str(MAOrd) '_FWHM0_' dtmethod{1} '_tVALUE_PW.nii.gz'];
%                     if ~exist(A,'file'); disp(['doesnot exists: ' A]); continue; end; 

                tstat          = alltstat{dt_cnt}{pwm_cnt}{AROcnt}{s_cnt};
                df             = T-rankXlist(dt_cnt);

                [p_ppx(s_cnt,:),p_ppy(s_cnt,:)]  = llpplot(tstat,df);
                [n_ppx(s_cnt,:),n_ppy(s_cnt,:)]  = llpplot(-tstat,df);

            end

            p_sappx(:,AROcnt) = mean(p_ppx);
            p_sappy(:,AROcnt) = mean(p_ppy);
            
            n_sappx(:,AROcnt) = mean(n_ppx);
            n_sappy(:,AROcnt) = mean(n_ppy);            
        end

        spfhRow1=subplot(2,3,dt_cnt);
        hold on; grid on; box on; 
        title([ 'T' num2str(T) ', ' pwmethod{1} ', ' dtmethod{1} ', Right Tail'],'Interpreter','latex')
        plot(-log10(p_sappx),-log10(p_sappy),'LineWidth',1.3)
        
        XA=get(gca,'Xlim');
        plot(XA(1):XA(2),XA(1):XA(2),'LineStyle','-.','color',[.6 .6 .6])
        
        line([cvt_05 cvt_05],XA,'LineWidth',1.1,'LineStyle','-.','color','r')
        line([cvt_01 cvt_01],XA,'LineWidth',1.1,'LineStyle','-.','color','r')
        line([cvt_001 cvt_001],XA,'LineWidth',1.1,'LineStyle','-.','color','r')
        
        xlabel('-log10(1/)','Interpreter','latex')
        ylabel('-log10(p-values)','Interpreter','latex')

        legend(cellfun(@num2str,num2cell(AROlist),'UniformOutput',false),'Location','northwest');
        
        spfhRow2=subplot(2,3,dt_cnt+3);
        hold on; grid on; box on; 
        title([ 'T' num2str(T) ', ' pwmethod{1} ', ' dtmethod{1} ', Left Tail'],'Interpreter','latex')
        plot(-log10(n_sappx),-log10(n_sappy),'LineWidth',1.3)
        
        XA=get(gca,'Xlim');
        plot(XA(1):XA(2),XA(1):XA(2),'LineStyle','-.','color',[.6 .6 .6])
        
        line([cvt_05 cvt_05],XA,'LineWidth',1.1,'LineStyle','-.','color','r')
        line([cvt_01 cvt_01],XA,'LineWidth',1.1,'LineStyle','-.','color','r')
        line([cvt_001 cvt_001],XA,'LineWidth',1.1,'LineStyle','-.','color','r')
        
        xlabel('-log10(1/)','Interpreter','latex')
        ylabel('-log10(p-values)','Interpreter','latex')

        legend(cellfun(@num2str,num2cell(AROlist),'UniformOutput',false),'Location','northwest');
        
        dt_cnt = dt_cnt + 1; 
    end

    %clear sadFPR
    pwm_cnt = pwm_cnt + 1; 
    set(psfh,'color','w')
    export_fig(psfh,[FigDir '/Fig_ppplots_T' num2str(T) '_' pwmethod{1} '_FWHM0.png'])

end
