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