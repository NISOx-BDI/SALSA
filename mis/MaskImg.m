function [S,MaskIdx] = MaskImg(Y,mask,ImgStat)
%[S,MaskIdx] = MaskImg(Y,mask,ImgStat)
%
% Y       : 2D Matrix OR path to the image. 
% mask    : the mask which should be applied on the brain 
% ImgStat : structure which comes from CleanNIFTI_spm.m
%
% S       : Cell , containing the time series for each seg/parcel
% MaskIdx : Cell , containing indices in the bigger 2D matrix
%
%
% Take the input as form of a matrix or string
% Extract the mask bit, calculates the mask and returns the masked matrix
%
% SA, Ox, 2020
verbose = 1; 
if ischar(Y) && nargin==2
    [Y,ImgStat] = CleanNIFTI_spm(Y,'verbose',0); 
    X0                          = ImgStat.ImgDim(1); 
    Y0                          = ImgStat.ImgDim(2); 
    Z0                          = ImgStat.ImgDim(3);
    I00                         = prod([X0,Y0,Z0]);
elseif isnumeric(Y) && nargin==3
    X0                          = ImgStat.ImgDim(1); 
    Y0                          = ImgStat.ImgDim(2); 
    Z0                          = ImgStat.ImgDim(3);
    T0                          = ImgStat.ImgDim(4);
    I00                         = prod([X0,Y0,Z0]);
    img_idx                     = 1:I00;
    img_idx(ImgStat.Removables) = []; % this is index of signal
    
    Y_tmp                       = zeros(I00,T0); %leave the T0 here, if it is one it is fine in case of 3D images
    Y_tmp(img_idx,:)            = Y;
    Y                           = reshape(Y_tmp,[X0 Y0 Z0 T0]); %leave the T0 here, if it is one it is fine in case of 3D images
else
    error('MaskImg:: Check inputs.')        
end

if ischar(mask)
    
    CLK	 = fix(clock);
    tmpdir  = [tempdir 'octspm12/tmp_' num2str(randi(5000)) '_' num2str(CLK(end))]; % make a temp directory 
    mkdir(tmpdir)
    if verbose; disp(['MaskImg:: --created: ' tmpdir]); end;
    
    if verbose; disp('MaskImg:: - gunzip the mask.'); end; 
    masktmp=[tmpdir '/mask_tmp_' num2str(randi(50)) num2str(CLK(end)+randi(1000)) '.nii'];
    system(['gunzip -c ' mask ' > ' masktmp]);    
    
    Vmask             = spm_vol(masktmp);
    mask              = spm_read_vols(Vmask);
    
    
    status = system(['rm -r ' tmpdir]);
    if ~status
        if verbose; disp(['--removed: ' tmpdir]); end; 
    else
        disp('MaskImg:: Warning: the temp directory was not deleted.')
    end    
    
else
    error('MaskImg:: mask should be a string; path 2 the mask')
end

umask                   = unique(mask);
umask(umask==0)         = []; 
if verbose; disp(['MaskImg:: The unique values in the mask: ' num2str(numel(umask))]); end; 


for m = umask'
    if verbose; disp(['MaskImg:: extracting: ' num2str(m)]); end; 
    tmask                    = mask; 
    tmask(tmask~=m)          = 0;
    tmask(tmask>0)           = 1; % binarise the mask temprorily

    tS                       = Y.*tmask; 
    tS                       = reshape(tS,I00,T0);
    tS(ImgStat.Removables,:) = []; 
    %sumS                    = round(sum(S,2),5); 
    %sumS                     = round(sum(tS,2).*10e5)./10e5; % Octave round function is dumb
    %tMaskIdx                 =  sumS>0; 
    tVar                     = var(tS,0,2);
    tMaskIdx                 =  tVar~=0; 
    tS(~tMaskIdx,:)          = [];
    
    S{m}                     = tS;
    MaskIdx{m}               = tMaskIdx;
end

end