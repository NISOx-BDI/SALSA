
CohortID = 'tHCP';
TaskName = 'GAMBLING';
T  = 253; 
TR = 0.72;

path2saveEVs = ['/Users/sorooshafyouni/Home/GitClone/FILM2/mis/EVs/' CohortID];

ParadNameList = {'loss_event','win_event','loss','win'};

figure; hold on; 
for tt = ParadNameList
    filename = ['/Users/sorooshafyouni/Home/GitClone/FILM2/mis/EVs/' CohortID '/' TaskName '/' tt{1} '.txt'];
    [EXD,EventTrail,OnSet,Duration,Events] = readHCPEvent0(filename,T,TR);
    
    
    filename=[path2saveEVs '/' TaskName '/task-' TaskName '_acq-' num2str(TR*1000) '_' tt{1} '_event.txt'];
    fileID = fopen(filename,'w');
    fprintf(fileID,'%f\n',EXD');
    fclose(fileID);    
    
    plot(EXD); 
    
end
legend(ParadNameList)

function [EDX,EventTrail,OnSet,Duration,Events] = readHCPEvent0(Path2Event,T,TR)

    dataArray = load(Path2Event);
    
    
    % Create output variable
    
    OnSet    = round(dataArray(:,1)./TR); 
    tOnSet   = OnSet;
    
    Duration = round(dataArray(:,2)./TR);  
    tDuration   = Duration;    
    
    EventTrail = zeros(round(T),1); 
    for i=1:numel(tOnSet)
        EventTrail(tOnSet(i):(tOnSet(i)+tDuration(i))) = 1;
    end
    
    hrf         = spm_hrf(TR);
    conv_dsgn1  = conv(EventTrail,hrf); 
    EDX         = conv_dsgn1(1:T);    
    
%     dataArray([1, 2]) = cellfun(@(x) num2cell(x), dataArray([1, 2]), 'UniformOutput', false);
%     Events = [dataArray{1:end-1}];
    Events = []; 
end