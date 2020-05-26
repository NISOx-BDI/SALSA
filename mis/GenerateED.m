
function EDX = GenerateED(lD,T,TR,onset)
% lD: length of stimuli in seconds
% T : length of the scan
% TR: TR in second

    tonset      = onset.*2; % use the same onset for the end too.
    Ted         = T-tonset;
    BCWidth     = fix(lD/2/TR);
    
    oneSTIM     = [zeros(1, BCWidth), ones(1, BCWidth)];
    
    numberSTIM  = fix(Ted/numel(oneSTIM));
    DD          = repmat(oneSTIM, [1, numberSTIM]);
    DD          = [zeros(1,onset) DD zeros(1,onset)];
    
    hrf         = spm_hrf(TR); 
    conv_dsgn   = conv(DD,hrf); 
    EDX         = conv_dsgn(1:T)';

end