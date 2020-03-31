function dY = detrend_polyfitM(Y,P,T)
% if numel(varargin)>=1
%     order = varargin{1};
% else
%     order = 2;
% end

sY = size(Y); 
if sY(2)~=T; Y = Y'; end
[nv,nt] = size(Y); 

for vi = 1:nv
    if ~mod(vi,10000); disp(['detrend_polyfit_voxel:: on voxel ' num2str(vi)]); end;         
    Ptmp = polyfitM((1:nt),Y(vi,:),P);
    dY(vi,:) = Y(vi,:)-polyval(Ptmp,(1:nt));
end


end





