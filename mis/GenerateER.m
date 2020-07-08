clear

CohortID           = 'ROCKLAND';
nsub               = 1;
NumTask            = 1; 

addpath('/Users/sorooshafyouni/Home/matlab/spm12') 

T                  = 404;
TR                 = 1.4;
BCl                = 10;

EDX = GenerateER0(T,TR,1,BCl,0); 

path2saveEVs = ['/Users/sorooshafyouni/Home/GitClone/FILM2/mis/EVs/' CohortID];
filename=[path2saveEVs '/' CohortID '_ERF_T' num2str(T) '_TR' num2str(TR*1000) '_taskpersec' num2str(BCl) '.txt'];
fileID = fopen(filename,'w');
fprintf(fileID,'%f\n',[EDX]');
fclose(fileID);


function EDX = GenerateER0(T,TR,nsub,TaskRate,FixedFlag)

    if ~exist('TaskRate','var'); TaskRate=5; end; 
        
    DinSec            = round(T*TR);
    DinTR             = T;
    hrf               = spm_hrf(TR);

    numEvents         = round(DinSec/TaskRate); % Keeping the rate as an event every 5seconds

    if ~FixedFlag
        minRest_inSec     = TaskRate/2; 
        maxRest_inSec     = TaskRate*2; 
    else
        minRest_inSec     = TaskRate;
        maxRest_inSec     = TaskRate;
    end
    
    minActivity_inSec = 1;
    maxActivity_inSec = 1;

    minRest_inTR      = fix(minRest_inSec/TR);
    maxRest_inTR      = fix(maxRest_inSec/TR);
    
    
    minActivity_inTR  = max(fix(minActivity_inSec/TR),1);
    maxActivity_inTR  = max(fix(maxActivity_inSec/TR),1);

    for sub_cnt = 1:nsub

        previousOnset    = 0;
        previousActivity = 0;
        temp1 = zeros(DinTR,1);
        
        for event = 1:numEvents
                
                if ~FixedFlag
                    EventRestinTR     = randi(maxRest_inTR);
                    EventActivityinTR = randi(maxActivity_inTR);
                else
                    EventRestinTR     = maxRest_inTR;
                    EventActivityinTR = maxActivity_inTR; 
                end
            
                randomOnsets     = previousOnset + previousActivity + max(EventRestinTR,minRest_inTR);
                randomDurations  = max(EventActivityinTR,minActivity_inTR);

                if randomOnsets > (T-(minRest_inTR*(3/TR))); continue; end; 

                
                for t = 0:(randomDurations-1)
                    temp1(randomOnsets+t) = 1;
                end

                previousOnset       = randomOnsets;
                previousActivity    = randomDurations;

        end
        conv_dsgn1     = conv(temp1,hrf); 
        EDX(:,sub_cnt) = conv_dsgn1(1:T)';
        
    end

end

% ff1=figure('position',[50,500,700,1000]);
% hold on; box on; axis tight; 
% title('Task2')
% imagesc(EDX)
% set(ff1,'color','w')    
% % export_fig('design_OneRegperEV.png')
% % 
% ff0=figure('position',[50,500,700,1000]);
% hold on; 
% for ii = 1:size(EDX,2)
%     subplot(size(EDX,2),1,ii); hold on; axis off
%     plot(EDX(:,ii))
% end
% set(ff0,'color','w')    
% export_fig('Reg_OneRegperEV.png')


