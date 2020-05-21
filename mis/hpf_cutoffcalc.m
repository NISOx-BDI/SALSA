function cutoff = hpf_cutoffcalc(X,TR,path2design)
% Calculates the cutoff for the highpass filter, given the deisgn
%
% SA,Ox,2020
%
    [ffpathstr,ffname,ffext] = fileparts(path2design);
    filename                 = [ffpathstr '/' ffname '.txt'];
    fileID                   = fopen(filename,'w');
    strpattern               = [repmat('%f',1,size(X,2)) '\n'];
    fprintf(fileID,strpattern,X);
    fclose(fileID);
    
    call_fsl(['${FSLDIR}/bin/Text2Vest ' filename ' ' path2design]);

    callcutoff          = ['${FSLDIR}/bin/cutoffcalc -i ' path2design ' --tr=' num2str(TR)];
    [status,cutoffstr]  = call_fsl(callcutoff);
    
    if ~status
        cutoff = str2double(cutoffstr);
    else
        error('the call to FSL was failed.')
    end
end