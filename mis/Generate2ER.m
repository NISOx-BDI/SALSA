function [EDX1,EDX2] = Generate2ER(T,TR,nsub,TaskRate)

    DinSec                  = round(T*TR);
    DinTR                   = T;
    hrf                     = spm_hrf(TR);
    
    numEvents               = round(DinSec/TaskRate);
    minRest_inSec           = 3; % 10
    maxRest_inSec           = 3; % 13

    minActivity_inSec       = 3; %3
    maxActivity_inSec       = 3; %7 

    minRest_inTR           = round(minRest_inSec/TR);
    maxRest_inTR           = round(maxRest_inSec/TR);
    minActivity_inTR       = round(minActivity_inSec/TR);
    maxActivity_inTR       = round(maxActivity_inSec/TR);

%     allRandomDurations1    = zeros(nsub,numEvents);
%     allRandomDurations2    = zeros(nsub,numEvents);
% 
%     allRandomOnsets1       = zeros(nsub,numEvents);
%     allRandomOnsets2       = zeros(nsub,numEvents);

    for sub_cnt = 1:nsub

        randomDurations1 = zeros(numEvents,1);
        randomDurations2 = zeros(numEvents,1);

        randomOnsets1    = zeros(numEvents,1);
        randomOnsets2    = zeros(numEvents,1);

        regressor1Events = 1;
        regressor2Events = 1;
        previousOnset    = 0;
        previousActivity = 0;

        for event = 1:(numEvents)

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
                randomOnsets1(regressor1Events)     = previousOnset + previousActivity + max(randi(maxRest_inTR),minRest_inTR);
                randomDurations1(regressor1Events)  = max(randi(maxActivity_inTR),minActivity_inTR);
                previousOnset       = randomOnsets1(regressor1Events);
                previousActivity    = randomDurations1(regressor1Events);
                regressor1Events    = regressor1Events + 1;
            elseif regressor == 2
                randomOnsets2(regressor2Events)     = previousOnset + previousActivity + max(randi(maxRest_inTR),minRest_inTR);
                randomDurations2(regressor2Events)  = max(randi(maxActivity_inTR),minActivity_inTR);
                previousOnset       = randomOnsets2(regressor2Events);
                previousActivity    = randomDurations2(regressor2Events);
                regressor2Events    = regressor2Events + 1;
            end
        
        end

        temp1 = zeros(DinTR,1);
        temp2 = zeros(DinTR,1);

        for event = 1:numEvents
            for t = 0:(randomDurations1(event)-1)
                RONS = randomOnsets1(event)+t;
                if RONS > (T-(minRest_inTR*(3/TR))); continue; end; 
                temp1(RONS) = 1;
            end
            for t = 0:(randomDurations2(event)-1)
                RONS = randomOnsets2(event)+t;
                if RONS > (T-(minRest_inTR*(3/TR))); continue; end; 
                temp2(RONS) = 1;
            end
        end
%       size(temp1)
%       figure; plot(temp1)
%       stims1(:,sub_cnt) = temp1;
%       stims2(:,sub_cnt) = temp2; 

        conv_dsgn1      = conv(temp1,hrf); 
        EDX1(:,sub_cnt) = conv_dsgn1(1:T)';

        conv_dsgn2      = conv(temp2,hrf); 
        EDX2(:,sub_cnt) = conv_dsgn2(1:T)';    

%       figure(1)
%       plot(temp1,'r')
%       hold on
%       plot(temp2,'b')
%       hold off
%       drawnow
%       pause(0.1)

%       allRandomOnsets1(sub_cnt,:) = randomOnsets1;
%       allRandomOnsets2(sub_cnt,:) = randomOnsets2;
% 
%       allRandomDurations1(sub_cnt,:) = randomDurations1;
%       allRandomDurations2(sub_cnt,:) = randomDurations2;

    end
end

