clear


addpath('/Users/sorooshafyouni/Home/matlab/spm12')

CohortID='ROCKLAND';

Path2Mats=['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_' CohortID '/']; 
FigDir = '/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/Img/FPR/RFigs';
SubIDPath = [Path2Mats '/' CohortID '_subid.mat'];
load(SubIDPath)
SubList=participants(1:51); 

nsub               = 50;
T                  = 900;
TR                 = 0.645;

DinSec            = round(T*TR);
DinTR             = T;
hrf               = spm_hrf(TR);

numEvents         = 18;
minRest_inSec     = 10;
maxRest_inSec     = 13;
minActivity_inSec = 3;
maxActivity_inSec = 7;

minRest           = fix(minRest_inSec/TR);
maxRest           = fix(maxRest_inSec/TR);
minActivity       = fix(minActivity_inSec/TR);
maxActivity       = fix(maxActivity_inSec/TR);


allRandomDurations1 = zeros(nsub,numEvents);
allRandomDurations2 = zeros(nsub,numEvents);

allRandomOnsets1 = zeros(nsub,numEvents);
allRandomOnsets2 = zeros(nsub,numEvents);

path2saveEVs='/Users/sorooshafyouni/Home/GitClone/FILM2/mis/EVs';

allsts = zeros(nsub,1);

for sub_cnt = 1:nsub
    
    randomDurations1 = zeros(numEvents,1);
    randomDurations2 = zeros(numEvents,1);
    
    randomOnsets1    = zeros(numEvents,1);
    randomOnsets2    = zeros(numEvents,1);
    
    regressor1Events = 1;
    regressor2Events = 1;
    previousOnset    = 0;
    previousActivity = 0;
    
    for event = 1:(numEvents*2)
        
        % All regressor 1 events have been set, the rest is regressor 2
        if (regressor1Events > numEvents)
            regressor = 2;
            % All regressor 2 events have been set, the rest is regressor 1
        elseif (regressor2Events > numEvents)
            regressor = 1;
        else
            % Randomly pick regressor 1 or 2
            regressor = randi(2);
        end
        
        if regressor == 1
            randomOnsets1(regressor1Events)     = previousOnset + previousActivity + max(randi(maxRest),minRest);
            randomDurations1(regressor1Events)  = max(randi(maxActivity),minActivity);
            previousOnset       = randomOnsets1(regressor1Events);
            previousActivity    = randomDurations1(regressor1Events);
            regressor1Events    = regressor1Events + 1;
        elseif regressor == 2
            randomOnsets2(regressor2Events)     = previousOnset + previousActivity + max(randi(maxRest),minRest);
            randomDurations2(regressor2Events)  = max(randi(maxActivity),minActivity);
            previousOnset       = randomOnsets2(regressor2Events);
            previousActivity    = randomDurations2(regressor2Events);
            regressor2Events    = regressor2Events + 1;
        end
        
    end
    
    sts(sub_cnt) = max(max(randomOnsets1(:)),max(randomOnsets2(:)));
    
    temp1 = zeros(DinTR,1);
    temp2 = zeros(DinTR,1);
    
    for event = 1:numEvents
        for t = 0:(randomDurations1(event)-1)
            temp1(randomOnsets1(event)+t) = 1;
        end
        for t = 0:(randomDurations2(event)-1)
            temp2(randomOnsets2(event)+t) = 1;
        end
    end
    
    stims1(:,sub_cnt) = temp1;
    stims2(:,sub_cnt) = temp2; 
    
    conv_dsgn1      = conv(temp1,hrf); 
    EDX1(:,sub_cnt) = conv_dsgn1(1:T)';
    
    conv_dsgn2      = conv(temp2,hrf); 
    EDX2(:,sub_cnt) = conv_dsgn2(1:T)';    
 
%     figure(1)
%     plot(temp1,'r')
%     hold on
%     plot(temp2,'b')
%     hold off
%     drawnow
%     pause(0.1)

    allRandomOnsets1(sub_cnt,:) = randomOnsets1;
    allRandomOnsets2(sub_cnt,:) = randomOnsets2;
    
    allRandomDurations1(sub_cnt,:) = randomDurations1;
    allRandomDurations2(sub_cnt,:) = randomDurations2;
    
%     filename=[path2saveEVs '/' CohortID '_sub_' SubList{sub_cnt} '_T' num2str(T) '_TR' num2str(TR*1000) '.txt'];
%     fileID = fopen(filename,'w');
%     fprintf(fileID,'%f %f\n',[EDX1(:,sub_cnt),EDX2(:,sub_cnt)]');
%     fclose(fileID);
    
end

psfh = figure('position',[50,500,650,600]);
subplot(3,1,1)
hold on; box on; grid on;
title('EV1')
plot(EDX1(:,1),'LineWidth',1.3)
plot(stims1(:,1),'linewidth',1.3,'linestyle','-.')
xlim([0 T])
legend({'HRF','Design'})
subplot(3,1,2)
hold on; box on; grid on; 
title('EV2')
plot(EDX2(:,1),'LineWidth',1.3)
plot(stims2(:,1),'linewidth',1.3,'linestyle','-.')
xlabel('time point')
xlim([0 T])
subplot(3,1,3)
hold on; grid on; box on; 
title('PSD(EV1-EV2)')
[xx,yy] = DrawMeSpectrum((EDX1(:,1)-EDX2(:,1)),TR);
plot(xx,yy,'LineWidth',1.3)
ylabel('Power')
xlabel('Frequency')

%set(psfh,'color','w')
%export_fig(psfh,['example_design.pdf'])
