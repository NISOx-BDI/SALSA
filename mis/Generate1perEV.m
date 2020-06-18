clear

nsub               = 100;

addpath('/Users/sorooshafyouni/Home/matlab/spm12') 
CohortID='ROCKLAND';

Path2Mats=['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_' CohortID '/']; 
FigDir = '/Users/sorooshafyouni/Home/GitClone/FILM2/NullRealfMRI/Img/FPR/RFigs';
SubIDPath = [Path2Mats '/' CohortID '_subid.mat'];
load(SubIDPath)
SubList = participants(1:nsub); 

% NKI 1400
% T                  = 404;
% TR                 = 1.4;

% NKI 900
% T                  = 900;
% TR                 = 0.645;

T                 = 900;
TR                = 0.645;

DinSec            = round(T*TR);
DinTR             = T;
hrf               = spm_hrf(TR);

numEvents         = 19;
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

path2saveEVs=['/Users/sorooshafyouni/Home/GitClone/FILM2/mis/EVs/RegpEV_' CohortID];
system(['mkdir -p ' path2saveEVs])

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
            randomOnsets1     = previousOnset + previousActivity + max(randi(maxRest),minRest);
            randomDurations1  = max(randi(maxActivity),minActivity);
            
            temp1 = zeros(DinTR,1);
            for t = 0:(randomDurations1-1)
                temp1(randomOnsets1+t) = 1;
            end
            
            conv_dsgn1               = conv(temp1,hrf); 
            EDX1(:,regressor1Events) = conv_dsgn1(1:T)';
            
            previousOnset       = randomOnsets1;
            previousActivity    = randomDurations1;
            regressor1Events    = regressor1Events + 1;
            
        elseif regressor == 2
            randomOnsets2     = previousOnset + previousActivity + max(randi(maxRest),minRest);
            randomDurations2  = max(randi(maxActivity),minActivity);
            
            temp2 = zeros(DinTR,1);
            for t = 0:(randomDurations2-1)
                temp2(randomOnsets2+t) = 1;
            end        
            
            conv_dsgn2               = conv(temp2,hrf); 
            EDX2(:,regressor2Events) = conv_dsgn2(1:T)'; 
            
            previousOnset       = randomOnsets2;
            previousActivity    = randomDurations2;
            regressor2Events    = regressor2Events + 1;
        end
                        
    end
    EDX = [EDX1,EDX2];
    filename=[path2saveEVs '/RegpEV_' CohortID '_sub_' SubList{sub_cnt} '_T' num2str(T) '_TR' num2str(TR*1000) '.txt'];
    fileID = fopen(filename,'w');
    strform = [repmat('%f ',1,size(EDX,2)) '\n'];
    fprintf(fileID,strform,[EDX1,EDX2]');
    fclose(fileID);
    
end


% ff1=figure('position',[50,500,700,1000]);
% subplot(1,2,1); hold on; box on; axis tight; 
% title('Task1')
% imagesc(EDX1)
% subplot(1,2,2); hold on; box on; axis tight; 
% title('Task2')
% imagesc(EDX2)
% set(ff1,'color','w')    
% export_fig('design_OneRegperEV.png')
% 
% EDX = [EDX1,EDX2];
% ff0=figure('position',[50,500,700,1000]);
% hold on; 
% for ii = 1:size(EDX,2)
%     subplot(size(EDX,2),1,ii); hold on; axis off
%     plot(EDX(:,ii))
% end
% set(ff0,'color','w')    
% export_fig('Reg_OneRegperEV.png')


