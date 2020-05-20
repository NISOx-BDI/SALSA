function ACL = ACLImage(Path2Img,Path2ACL)

    [Y,ImgStat]   = CleanNIFTI_spm(Path2Img);
    T = size(Y,2);
    acf = AC_fft(Y,T);
    ACL = sum(acf(:,1:round(T/4)).^2,2);
    CleanNIFTI_spm(ACL,'ImgInfo',ImgStat.spmV,'destdir',Path2ACL,'Removables',ImgStat.Removables);

end