clear

addpath('~/bin/FILM2/mis')

pwmethodlist={'ACF','AR-YW','AR-W','ARMAHR'};
AROlist={{30,60};{1,2,5,10,20};{1,2,5,10,20};{1,2,5,10,20}};
dtrendmethodlist={'spline','poly','dct','hpf'};
dtrendOlist=[NaN,4,10,100];
MAOlist = [0,0,0,1];
FWHMl=0;
subidlist={'A00008326','A00008399','A00010893','A00013809','A00018030','A00019903','A00021039','A00023510','A00025566','A00027159','A00027167','A00027439','A00027443','A00027544','A00027651','A00028150','A00028152','A00028177','A00028184','A00028185'};
sesid='DS2';
get_tcrit = @(alph,df) tinv(1-alph,df);

T           = 900;
rankXlist   = [31,30,36,26];
alphlev     = 0.05;

ResulDir='/well/nichols/users/scf915/ROCKLAND/R.PW';

pwm_cnt = 1; 
for pwmethod=pwmethodlist
    
    dtm_cnt = 1; 
    for dtm=dtrendmethodlist
       disp([pwmethod{1} ' - ' dtm{1} ' - ' num2str(dtrendOlist(dtm_cnt))])
       
       ar_cnt=1;
       for ARO=AROlist{pwm_cnt}
           ResultSubDir=[ResulDir '/ROCKLAND_' pwmethod{1} '_AR-' num2str(ARO{1}) '_MA-' num2str(MAOlist(pwm_cnt)) '_FWHM' num2str(FWHMl) '_' dtm{1}];
           
           s_cnt = 1;
           for subid=subidlist
             ResultSubjectDir=[ResultSubDir '/' subid{1} '_' sesid];
             if isnan(dtrendOlist(dtm_cnt)); dtmO=''; else; dtmO=num2str(dtrendOlist(dtm_cnt)); end
             ImageFileName=[ResultSubjectDir '/EDboxcar_20_' pwmethod{1} '_AR' num2str(ARO{1}) '_MA' num2str(MAOlist(pwm_cnt)) '_FWHM' num2str(FWHMl) '_' dtm{1} dtmO '_tVALUE_PW.nii'];
             if exist(ImageFileName,'file')
                 Y = CleanNIFTI_spm(ImageFileName,'verbose',0);

                 
                 
                 df      = T-rankXlist(dtm_cnt);
                 t_crit  = get_tcrit(alphlev/2,df);
                 
                 num_sigvox  = sum(abs(Y)>t_crit);
                 FPR(s_cnt,pwm_cnt,ar_cnt,dtm_cnt)         = num_sigvox./numel(Y);
                 s_cnt = s_cnt + 1;
             end
           end
           ar_cnt = ar_cnt + 1;
       end
        dtm_cnt = dtm_cnt + 1; 
    end
    pwm_cnt = pwm_cnt + 1; 
end


pwm_cnt = 1; 
for pwmethod=pwmethodlist
    dtm_cnt = 1; 
    for dtm=dtrendmethodlist
       ar_cnt=1;
       for ARO=AROlist{pwm_cnt}           
                 sFPR(pwm_cnt,ar_cnt,dtm_cnt)=mean(nonzeros(FPR(:,pwm_cnt,ar_cnt,dtm_cnt)));
           ar_cnt = ar_cnt + 1;
       end
        dtm_cnt = dtm_cnt + 1; 
    end
    pwm_cnt = pwm_cnt + 1; 
end

ARORDER={1,2,5,10,20};
fh = figure; 
for i = 1:5
    sp=subplot(3,2,i);
    hold on; box on; grid on; 
    title(['Model Order:' num2str(ARORDER{i})])
    bar(squeeze(sFPR(:,i,:)))
    ylim([0 0.1])
    plot(0:5,0.05.*ones(1,6))
    
    sp.XTick=1:4;
    sp.XTickLabel=pwmethodlist;
    
    legend(dtrendmethodlist,'location','northwest')
    
end
    
set(fh,'color','w')


