function ZStat = myt2z(Path2T,Path2Z,df)

    [TStat,ImgStat]   = CleanNIFTI_spm(Path2T,'verbose',0);
    ZStat    = spm_t2z(TStat,df);
    CleanNIFTI_spm(ZStat,'ImgInfo',ImgStat.spmV,'destdir',Path2Z,'Removables',ImgStat.Removables,'verbose',0);

end