% clear
% 
% CohortID           = 'ROCKLAND';
% nsub               = 1;
% NumTask            = 1; 
% 
% addpath('/Users/sorooshafyouni/Home/matlab/spm12') 
% 
% T                  = 900;
% TR                 = 0.645;
% 
% EDX = GenerateORPE(T,TR,1); 

function EDX = GenerateER(T,TR,nsub,TaskRate)

    if ~exist('TaskRate','var'); TaskRate=5; end; 
        
    DinSec            = round(T*TR);
    DinTR             = T;
    hrf               = spm_hrf(TR);

    numEvents         = round(DinSec/TaskRate); % Keeping the rate as an event every 5seconds

    minRest_inSec     = TaskRate/2; 
    maxRest_inSec     = TaskRate*2; 

    minActivity_inSec = 1;
    maxActivity_inSec = 1;

    minRest_inTR      = fix(minRest_inSec/TR);
    maxRest_inTR      = fix(maxRest_inSec/TR);
    minActivity_inTR  = fix(minActivity_inSec/TR);
    maxActivity_inTR  = fix(maxActivity_inSec/TR);

    for sub_cnt = 1:nsub

        previousOnset    = 0;
        previousActivity = 0;
        temp1 = zeros(DinTR,1);
        
        for event = 1:(numEvents)

                randomOnsets     = previousOnset + previousActivity + max(randi(maxRest_inTR),minRest_inTR);
                randomDurations  = max(randi(maxActivity_inTR),minActivity_inTR);

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


