function filter=hp_fsl(Nscans,cut,TR)
%This function creates the highpass filter matrix that is used in FSL for a
%given highpass filter cutoff.  To apply the filter to a timeseries simply
%multiply the filter and the data (filter goes on the left side of the
%multiplication).
%
%Nscans:  The number of time points in the time series you are filtering
%cut:  The highpass filter cutoff in seconds
%TR:  The TR for your data in seconds


cut=cut/TR;
 sigN2=(cut/sqrt(2))^2;
 K    = toeplitz(1/sqrt(2*pi*sigN2)*exp(-[0:(Nscans - 1)].^2/(2*sigN2)));
 K    = spdiags(1./sum(K')',0,Nscans,Nscans)*K;
 
 H = zeros(Nscans,Nscans); % Smoothing matrix, s.t. H*y is smooth line
 X = [ones(Nscans,1) (1:Nscans)'];
 for k = 1:Nscans
   W = diag(K(k,:));
   Hat = X*pinv(W*X)*W;
   H(k,:) = Hat(k,:);
 end
 
 filter=eye(Nscans)-H;