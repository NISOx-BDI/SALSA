% filename = '/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/ROCKLAND/R_mpp/sub-A00028150/ses-DS2/sub-A00008399_ses-DS2_task-CHECKERBOARD_acq-645_bold_mpp/sub-A00008399_ses-DS2_task-CHECKERBOARD_acq-645_events.tsv';
% 
% 
% T  = 260; 
% TR = 0.645;
% [EXD,EventTrail,OnSet,Duration,Events] = readBIDSEvent0(filename,T,TR);

function [EDX,EventTrail,OnSet,Duration,Events] = readBIDSEvent(Path2Event,T,TR)

    delimiter = '\t';
    startRow = 2;

    % Format for each line of text:
    %   column1: double (%f)
    %	column2: double (%f)
    %   column3: text (%s)
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%f%f%s%[^\n\r]';

    % Open the text file.
    fileID = fopen(Path2Event,'r');

    % Read columns of data according to the format.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    
    % Close the text file.
    fclose(fileID);

    % Post processing for unimportable data.
    % No unimportable data rules were applied during the import, so no post
    % processing code is included. To generate code which works for
    % unimportable data, select unimportable cells in a file and regenerate the
    % script.

    
    % Create output variable
    
    OnSet    = round(dataArray{1}./TR); 
    rOnSet   = OnSet(1:2:end);
    tOnSet   = OnSet(2:2:end);
    
    Duration = round(dataArray{2}./TR);  
    sDuration   = Duration(1:2:end);
    tDuration   = Duration(2:2:end);    
    
    EventTrail = zeros(round(T),1); 
    for i=1:numel(tOnSet)
        EventTrail(tOnSet(i):(tOnSet(i)+tDuration(i))) = 1;
    end
    
    hrf         = spm_hrf(TR);
    conv_dsgn1  = conv(EventTrail,hrf); 
    EDX         = conv_dsgn1(1:T);    
    
    dataArray([1, 2]) = cellfun(@(x) num2cell(x), dataArray([1, 2]), 'UniformOutput', false);
    Events = [dataArray{1:end-1}];
    
end