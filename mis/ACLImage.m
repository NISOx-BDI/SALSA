function ACL = ACLImage(V,Path2ACL)
% V: should be either the path to an image or a 2-D matrix
%    if 2-D matrix, should have form of VxT

if nargin==2
    [Y,ImgStat]   = CleanNIFTI_spm(V);
    T = size(Y,2);
    acf = AC_fft(Y,T);
    ACL = sum(acf(:,1:round(T/4)).^2,2);
    CleanNIFTI_spm(ACL,'ImgInfo',ImgStat.spmV,'destdir',Path2ACL,'Removables',ImgStat.Removables);

elseif nargin==1 && isnumeric(V)
        T   = size(V,2); 
        acf = AC_fft(V,T);
        ACL = sum(acf(:,1:round(T/4)).^2,2);
end