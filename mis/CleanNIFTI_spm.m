function [Y,Stat]=CleanNIFTI_spm(V0,varargin)
%[Y,Stat]=CleanNIFTI_fsl(V0,varargin)
% Nomalise, demean and clean (remove zeros and NaNs) from fMRI data. 
%
% Output Y is a IxT matrix
% Trigger 'DestDir' with the destination to save the NIFTI file. 
%
% Main difference from CleanNIFTI is that this function uses FSL and is
% operational on Octave.
%_________________________________________________________________________
% Soroosh Afyouni, NISOx.org, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________


Steps=[]; verbose=1; DestDir=[]; md=[]; scl=[]; ImgStat=[]; fwhm = 0;
Path2Mask = []; 
datatype  = 'f';
%voxelsize = [2 2 2 1];
SaveFlag     = 0;

if sum(strcmpi(varargin,'verbose'))
   verbose      =   varargin{find(strcmpi(varargin,'verbose'))+1}; 
end
%
if sum(strcmpi(varargin,'removables'))
   Removables      =   varargin{find(strcmpi(varargin,'removables'))+1}; 
end
%
if sum(strcmpi(varargin,'destdir'))
   DestDir      =   varargin{find(strcmpi(varargin,'destdir'))+1};
end

if sum(strcmpi(varargin,'mask'))
   Path2Mask      =   varargin{find(strcmpi(varargin,'mask'))+1};
end

if sum(strcmpi(varargin,'fwhm'))
   fwhm      =   varargin{find(strcmpi(varargin,'fwhm'))+1};
end

if sum(strcmpi(varargin,'ImgInfo'))
   ImgStat      =   varargin{find(strcmpi(varargin,'ImgInfo'))+1};
   SaveFlag     = 1;
end

%if sum(strcmpi(varargin,'voxelsize'))
%   voxelsize      =   varargin{find(strcmpi(varargin,'voxelsize'))+1};
%end
%if sum(strcmpi(varargin,'ImgDim'))
%   ImgDim      =   varargin{find(strcmpi(varargin,'ImgDim'))+1};
%end
%if sum(strcmpi(varargin,'Removables'))
%   img_Removables      =   varargin{find(strcmpi(varargin,'Removables'))+1};
%end

% untested, so I leave it as it is for now
% if sum(strcmpi(varargin,'bptf'))
%     bp      =   varargin{find(strcmpi(varargin,'bptf'))+1};
%     if numel(bp<3); error('should be: [lowerbound upperbound TR]'); end;
% else
%     bp=[];
% end
%
if sum(strcmpi(varargin,'norm'))
   scl          =   varargin{find(strcmpi(varargin,'norm'))+1};
end
%
if sum(strcmpi(varargin,'scale'))
   scl          =   varargin{find(strcmpi(varargin,'scale'))+1};
   md           =   1;
end

%
if sum(strcmpi(varargin,'demean'))
    dmflag          =   1;
else
    if verbose; warning('the demean flag is off!'); end;
    dmflag          =   0;
end

%sort out the FSL directories ==============================
% fsldir = getenv('FSLDIR');
% if isempty(fsldir) 
%     error('CleanNIFTI_fsl:: I can not find FSL directory!'); 
% end
% fsldirmpath = sprintf('%s/etc/matlab',fsldir);
% path(path, fsldirmpath);
% clear fsldir fsldirmpath;

CLK	 = fix(clock);
tmpdir  = [tempdir 'octspm12/tmp_' num2str(randi(5000)) '_' num2str(CLK(end))]; % make a temp directory 
mkdir(tmpdir)
if verbose; disp(['CleanNIFTI_spm:: --created: ' tmpdir]); end;

if ischar(V0)
        
    [ffpathstr,ffname,ffext]=fileparts(V0);
    if verbose; disp(['CleanNIFTI_spm:: -Path to the image is: ' ffpathstr]); end;
    
    if ~isempty(strfind(ffname,'.dtseries')) || ~isempty(strfind(ffext,'.dtseries'))
        if verbose; disp(['CleanNIFTI_spm:: --File is CIFTI: ' ffname ffext]); end;
        error('CleanNIFTI_spm:: CIFTI not supported yet.')
%        error('We can not support CIFTI files yet.')
%         V1=ft_read_cifti(V0);
%         V2=V1.dtseries;
%         I0=size(V2,1); T0=size(V2,2);
%         Y=V2; clear V2 V1; 
    elseif isempty(strfind(ffname,'.dtseries')) || ~isempty(strfind(ffname,'.nii'))
        if verbose; disp(['CleanNIFTI_spm:: --File is NIFTI: ' ffname ffext]); end;
        
        if strfind(V0,'.gz')
            if verbose; disp('CleanNIFTI_spm:: - gunzip the image.'); end; 
            randtempfilename=[tmpdir '/img_tmp_' num2str(randi(50)) num2str(CLK(end)+randi(1000)) '.nii'];
            system(['gunzip -c ' V0 ' > ' randtempfilename]);        
        else
            randtempfilename = V0;
        end
        
        Vmask   = spm_vol(randtempfilename);
        V2        = spm_read_vols(Vmask);
        dims      = Vmask(1).private.dat.dim;
        TR        = Vmask(1).private.timing.tspace;
        voxelsize = [abs(diag(Vmask(1).mat(1:3,1:3)))' TR];
        Stat.spmV = Vmask;
        %[V2,dims,voxelsize] = read_avw(V0);
        
        ND=ndims(V2);
        
        if ND==3
            T0 = 1;
        elseif ND==4
            T0 = size(V2,4);
        else
            error('CleanNIFTI_spm:: The input should be either 3D or 4D.');
        end
        X0 = size(V2,1); Y0 = size(V2,2); Z0 = size(V2,3); 
        
        if sum(ismember(dims,[X0,Y0,Z0,T0]))~=ND
            error('CleanNIFTI_spm:: something is wrong with the dimensions!'); 
        end
        
        %figure; imagesc(squeeze(mean(V2(:,25,:,:),4)))
        
        % Guassian Kernel Smoothing 
        if fwhm
            disp(['CleanNIFTI_spm:: -Smooth the image with fwhm: ' num2str(fwhm)])
            disp(['CleanNIFTI_spm:: -- voxel sizes: ' num2str(voxelsize(1)) 'x' num2str(voxelsize(2)) 'x' num2str(voxelsize(3))])
            V2 = GaussianSmooth(V2,fwhm,round(voxelsize));
        end
        
        %figure; imagesc(squeeze(mean(V2(:,25,:,:),4)))
        
        I0 = prod([X0,Y0,Z0]);
        Y  = reshape(V2,[I0,T0]); clear V2;
    else
        error('CleanNIFTI_spm:: Unknown input image.')
    end
    
    if verbose; disp('CleanNIFTI_spm:: -Image loaded.'); end;
    
elseif isnumeric(V0)
    if numel(size(V0))==2 && ismember(1,size(V0))
        if verbose; disp('CleanNIFTI_spm:: -Input is a 1D Matrix.'); end;
        Y = double(V0);  
        I0= numel(Y);
        T0=1;
        ImageType = 1; 
    elseif numel(size(V0))==2 && ~ismember(1,size(V0))
        if verbose; disp('CleanNIFTI_spm:: -Input is a 2D Matrix.'); end;
        if  size(V0,1)<size(V0,2); error('CleanNIFTI_spm:: Matrix should be voxels X time: Transpose the input'); end;
        Y = double(V0);  
        I0= size(Y,1); T0 = size(Y,2);
        ImageType = 2;
    elseif numel(size(V0))==3
         if verbose; disp('CleanNIFTI_spm:: -Input is a 2D Matrix.'); end;
         error('CleanNIFTI_spm:: Has not developed anything for this yet!')
         ImageType = 3;
    elseif numel(size(V0))==4
        disp(['CleanNIFTI_spm:: -Input is a 4D matrix.'])
        X0 = size(V0,1); Y0 = size(V0,2); Z0 = size(V0,3); T0 = size(V0,4);
        I0 = prod([X0,Y0,Z0]);
        Y  = reshape(V0,[I0,T0]);
        ImageType = 4;
    else
        error('CleanNIFTI_spm:: Something does not match with input dimensions.')
    end
end

%Y = double(Y);%to work with int 16bit as well.
if ~SaveFlag
    
    if ~isempty(Path2Mask)
        if verbose; disp(['CleanNIFTI_spm:: --A mask is being used to extract ''good'' time series.']); end; 
        if strfind(Path2Mask,'.gz')
            disp('CleanNIFTI_spm:: - gunzip the mask.')
            randtempfilename4mask=[tmpdir '/mask_tmp_' num2str(randi(50)) num2str(CLK(end)+randi(1000)) '.nii'];
            system(['gunzip -c ' Path2Mask ' > ' randtempfilename4mask]);
        else
            randtempfilename4mask = Path2Mask;
        end
        Vmask   = spm_vol(randtempfilename4mask);
        mask    = spm_read_vols(Vmask);
        maskdim = size(mask);
        mask    = reshape(mask,[1,prod(maskdim)]);
        Removables = find(mask==0);
    elseif isempty(Path2Mask)
        if verbose; disp(['-Any time series which is either nan or variance zero is removed.']); end
        nan_idx    = find(isnan(sum(Y,2)));
        %zeros_idx  = find(sum(abs(Y),2) < (eps*10e5) ); %find(sum(Y,2)==0);
        zeros_idx  = find(var(Y,0,2)==0);
        Removables = [nan_idx;zeros_idx];
    end
        
    idx        = 1:I0;
    idx(Removables) = [];
    Y(Removables,:) = [];
    I1 = size(Y,1); %update number of voxels

    if verbose; disp(['CleanNIFTI_spm:: -Extra-cranial areas removed: ' num2str(size(Y,1)) 'x' num2str(size(Y,2))]); end;
    Steps=[Steps 'CLEANED_'];

    Stat.GlobalMeanSignal = mean(Y);
    Stat.OrigDim     = [I0 T0];
    Stat.CleanedDim  = [I1 T0];
    Stat.Removables  = Removables;
    Stat.idx         = idx;
    OrigMean         = mean(Y,2);
    Stat.OrigMean    = OrigMean;
    Stat.ImgDim      = [X0 Y0 Z0 T0];
    [X3d,Y3d,Z3d]    = ind2sub(Stat.ImgDim(1:3),idx);
    Stat.idx3D       = [X3d;Y3d;Z3d]';       
    Stat.voxelsize   = voxelsize;
    Stat.datatype    = datatype;
    Stat.TR          = TR;
end
%------------------------------------------------------------------------
% Intensity Normalisation------------------------------------------------------
IntnstyScl = @(Y,md,scl) (Y./md).*scl; 
if ~isempty(scl) && isempty(md) && ~SaveFlag
    md  = median(mean(Y,2)); %NB median of the mean image.
    %md  = mean(mean(Y,2)); %NB *mean* of the mean image.
    Y   = IntnstyScl(Y,md,scl);
    if verbose; disp(['-Intensity Normalised by ' num2str(scl) '&' num2str(md) '.']); end;
    Steps = [Steps 'NORM_'];
elseif ~isempty(scl) && ~isempty(md)
    assert(md==1,'4D mean in scalling cannot be anything other than 1!')
    Y   = IntnstyScl(Y,md,scl);
    if verbose; disp(['-Intensity Scaled by ' num2str(scl) '.']); end;
    Steps = [Steps 'SCALE_'];
elseif isempty(scl) && isempty(md)    
    if verbose; disp('-No normalisation/scaling has been set!'); end;
else
    error('CleanNIFTI_spm:: Something is wrong with param re: intensity normalisation')
end
%------------------------------------------------------------------------
%Centre the data-----------------------------
if dmflag && ~SaveFlag
    mvY_NormInt      = mean(Y,2); %later will be used as grand mean! don't touch it!
    dmeaner          = repmat(mvY_NormInt,[1,T0]);
    Y                = Y-dmeaner; clear dmeaner
    DemeanedMen      = mean(Y,2);
    Stat.DemeanedMen = DemeanedMen;
    Steps            = [Steps 'DEMEANED_'];
end
%------------------------------------------------------------------------
%Save the image-----------------------------
if ~isempty(ImgStat) && isnumeric(V0) && SaveFlag
    %if ~any(strfind(path,'spm')); warning('**SPM has not been added to the path!**'); end;
    %DestDir = ImgStat(1).fname; 
    [spathstr,sname,stext]=fileparts(DestDir);
    
    
    if isempty(spathstr)
        Dir2Save = sname;
    else
        if exist(spathstr,'dir')~=7; mkdir(spathstr); end;
        Dir2Save = [spathstr '/' sname stext];
    end
    
    if verbose; disp(['CleanNIFTI_spm:: Image saved: ' Dir2Save]); end; 
    
    X0= ImgStat(1).dim(1); 
    Y0= ImgStat(1).dim(2); 
    Z0= ImgStat(1).dim(3);
    
    I00 = prod([X0,Y0,Z0]);
    
    if T0==1 
        disp('CleanNIFTI_spm:: Saved image will be 3D.');
        ImgStat = ImgStat(1); 
    end
    
    img_idx = 1:I00;
    img_idx(Removables) = []; % this is index of signal
    
    Y_tmp            = zeros(I00,T0); %leave the T0 here, if it is one it is fine in case of 3D images
    Y_tmp(img_idx,:) = Y;
    Y_tmp            = reshape(Y_tmp,[X0 Y0 Z0 T0]); %leave the T0 here, if it is one it is fine in case of 3D images
    
    for it = 1:T0
        %ImgStat(it).dim   = size(Y_tmp); 
        ImgStat(it).fname = Dir2Save;
        ImgStat(it).private.dat.fname = Dir2Save;
        spm_write_vol(ImgStat(it),Y_tmp(:,:,:,it));
    end
    sstmp = system(['gzip -f ' Dir2Save]);
    if sstmp 
        error('CleanNIFTI_spm:: Failed to save the file'); 
    end
    %save_avw(Y_tmp,Dir2Save,datatype,ImgStat.voxelsize);
    clear *_tmp clear V_Img;
else
    if verbose
        disp('-the images will NOT be saved:')
        disp('-- Either destination directory was not set OR the input is not a nifti.')
    end
    clear V_Img;
end

%------------------------------------------------------------------------ 
%Temporal filtering-----------------------------
% if ~isempty(bp)
%     disp(['bandpass filter of lowerbound ' num2str(bp(1)) ' & higherbound ' um2str(bp(2)) ' is runing!'])
%     fsl_bptf(DestDir,DestDir,[bp(1) bp(2),TR])
%     Steps=[Steps 'BPTF'];
% end
Stat.Steps = Steps;


status = system(['rm -r ' tmpdir]);
if ~status
    if verbose; disp(['--created: ' tmpdir]); end; 
else
    disp('CleanNIFTI_spm:: Warning: the temp directory was not deleted.')
end

end

% function fsl_bptf(Vin,Vout,bp,TR)
% % bp should be in frq
% f2tr=@(f,TR) 0.5*(1/f)/TR;
% system(['/usr/local/fsl/bin/fslmaths ' Vin ' -bptf ' num2str(round(f2tr(bp(1),TR),2)) ' ' num2str(round(f2tr(bp(2),TR),2)) ' ' Vout ])


function sV = GaussianSmooth(V,fwhm,D)
   % from FMRISTAT, gauss_blur.m

   numxs     = size(V,1);
   numys     = size(V,2);
   nZ        = size(V,3);
   nT        = size(V,4);
   numpix    = numxs*numys;
   
   if length(fwhm) == 1
      fwhm=repmat(fwhm,1,3);
   end
   fwhm_x = fwhm(1)/abs(D(1));
   ker_x  = exp(-(-ceil(fwhm_x):ceil(fwhm_x)).^2*4*log(2)/fwhm_x^2);
   ker_x  = ker_x/sum(ker_x);
   fwhm_y = fwhm(2)/abs(D(2));
   ker_y  = exp(-(-ceil(fwhm_y):ceil(fwhm_y)).^2*4*log(2)/fwhm_y^2);
   ker_y  = ker_y/sum(ker_y);
   fwhm_z = fwhm(3)/abs(D(3));
   ker_z  = exp(-(0:(nZ-1)).^2*4*log(2)/fwhm_z^2);
   K      = toeplitz(ker_z);
   K      = K./(ones(nZ)*K);
   
   sV=zeros(numxs,numys,nZ,nT);
   for iT = 1:nT
      for iZ = 1:nZ
         temp(:,iZ) = reshape(conv2(ker_x,ker_y,V(:,:,iZ,iT),'same'),numpix,1);   
      end
      sV(:,:,:,iT)  = reshape(temp*K,[numxs numys nZ 1]);
   end
end