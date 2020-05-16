
function EDX = GenerateED(lD,T,TR)
% lD: length of stimuli in seconds
% T : length of the scan
% TR: TR in second

    BCWidth     = fix(lD/2/TR);
    oneSTIM     = [zeros(1, BCWidth), ones(1, BCWidth)];
    
    numberSTIM  = fix(T/numel(oneSTIM));
    DD          = repmat(oneSTIM, [1, numberSTIM]);
    
    
    hrf         = spm_hrf(TR); 
    conv_dsgn   = conv(DD,hrf); 
    EDX         = conv_dsgn(1:T)';

end