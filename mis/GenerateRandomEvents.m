clear


addpath('/Users/sorooshafyouni/Home/matlab/spm12')

% CohortID           = 'HCP';
% T                  = 1200;
% TR                 = 0.72;

CohortID           = 'NEO';
T                  = 2300;
TR                 = 0.392;

nsub               = 51;

% NKI 1400
% T                  = 404;
% TR                 = 1.4;

% Beijing 225
% T                  = 225;
% TR               = 2;

TaskRate = 5;
[EDX1,EDX2] = Generate2ER(T,TR,nsub,TaskRate);

path2saveEVs=['/Users/sorooshafyouni/Home/GitClone/FILM2/mis/EVs/' CohortID];
system(['mkdir -p ' path2saveEVs])
Path2Mats=['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_' CohortID '/']; 
FigDir = '/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/Img/FPR/RFigs';
SubIDPath = [Path2Mats '/' CohortID '_subid.mat'];
load(SubIDPath)
SubList=participants(1:nsub); 
for sub_cnt = 1:nsub
    filename=[path2saveEVs '/' CohortID '_sub_' SubList{sub_cnt} '_T' num2str(T) '_TR' num2str(TR*1000) '.txt'];
    fileID = fopen(filename,'w');
    fprintf(fileID,'%f %f\n',[EDX1(:,sub_cnt),EDX2(:,sub_cnt)]');
    fclose(fileID);
end

psfh = figure('position',[50,500,650,600]);
subplot(3,1,1)
hold on; box on; grid on;
title('EV1')
plot(EDX1(:,1),'LineWidth',1.3)
%plot(stims1(:,1),'linewidth',1.3,'linestyle','-.')
xlim([0 T])
legend({'HRF'})
subplot(3,1,2)
hold on; box on; grid on; 
title('EV2')
plot(EDX2(:,1),'LineWidth',1.3)
%plot(stims2(:,1),'linewidth',1.3,'linestyle','-.')
xlabel('time point')
xlim([0 T])
subplot(3,1,3)
hold on; grid on; box on; 
title('PSD(EV1-EV2)')
[xx,yy] = DrawMeSpectrum((EDX1(:,1)-EDX2(:,1)),TR);
plot(xx,yy,'LineWidth',1.3)
ylabel('Power')
xlabel('Frequency')
figure; plot(EDX1); hold on; plot(EDX2)

%set(psfh,'color','w')
%export_fig(psfh,[CohortID '_T' num2str(T) '_TR' num2str(TR) '_example_design.pdf'])


