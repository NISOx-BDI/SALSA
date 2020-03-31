function [dymat,polycoefs,trendY] = multpolyfit(xmat,ymat,T,p)

% [dymat,polycoefs,trendY] = multpolyfit(xmat,ymat,p)
%
% 1) finds the coefficients of a polynomial P(X) of
%    degree p that fits the data Y best in a least-squares sense.
% 2) returns the value of a polynomial P evaluated at xmat
% 3) deternd ymat
%
% This is to fasten up the detrending. With Matlab function is takes couple
% of mins to finish on a 50K voxel image. 
% This should fasten up by 10x. 
% See the bottom for the Matlab implementation and the sanity check 
% 
% INPUTS:
% xmat: 
% ymat: time series in VxT -- V: number of voxels, T: length of time series
% p : order of the polynomial fit (well, I want to keep this function generic)
%
% OUTPUTS:
% dymat : detrended time series
% polycoeffs : polynomial coefficients in
%    ascending powers, FLIP[P(1)*X^N + P(2)*X^(N-1) +...+ P(N)*X + P(N+1)]
% trendY : Trend of the time series
%
% Soroosh Afyouni, Uni of Oxford, 2020


%%% FUNCTIONS
% substitute for polyval.m in Matlab -- ~5x faster
multpolyval = @(p, x) sum(shiftdim(p(:), -ndims(x)) .* x .^ shiftdim((numel(p)-1:-1:0).', -ndims(x)), ndims(x)+1);

%%% SCRIPT
% check stuff:
if any(~(size(xmat) == size(ymat)))
    error('xmat and ymat should have identical dimensions')
end
    
if  size(ymat,1) ~= T
    ymat = ymat'; xmat = xmat';
end

if nargin < 3 
    error('not enough input!')
elseif nargin < 4
   p = 2; % I wanted to implement AFNI 3dConvolve rule of thump here, but then I realised I need the TR... So will do it outside to keep the function generic for time series analysis  
end



[nl,nv]     = size(ymat);
disp(['multpolyfit:: image dims; ' num2str(nl) ' time x ' num2str(nv) ' voxels, Poly order: ' num2str(p)])

% to fasten up stuff
polycoefs   = zeros(nv,p+1);
trendY      = zeros(nv,nl); 
dymat       = zeros(nv,nl);

for iv = 1:nv
    if ~mod(iv,10000); disp(['multpolyfit:: on voxel: ' num2str(iv)]); end; 
    M       = repmat(xmat(:,iv),1,p+1);
    M       = bsxfun(@power,M,0:p);
    ctmp    = M\ymat(:,iv);
    polycoefs(iv,:) = ctmp; 
    tytmp   = multpolyval(flip(ctmp),xmat(:,iv));
    trendY(iv,:) = tytmp;  
    dymat(iv,:)  = ymat(:,iv)-tytmp;
end

end


%%% SANITY CHECK:
% c1 = zeros(n+1,m);
% for k=1:m
%    c1(:,k) = polyfit(x(:,k),y(:,k),n)';
% end
% This is what I did before, works perfectly, but it is really slow!
% sY = size(Y); 
% if sY(2)~=T; Y = Y'; end
% [nv,nt] = size(Y); 
% 
% for vi = 1:nv
%     if ~mod(vi,10000); disp(['detrend_polyfit_voxel:: on voxel ' num2str(vi)]); end;         
%     Ptmp = polyfit((1:nt),Y(vi,:),P);
%     dY(vi,:) = Y(vi,:)-polyval(Ptmp,(1:nt));
% end
