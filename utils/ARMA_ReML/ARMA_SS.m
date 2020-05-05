function [a,b,LL] = ARMA_SS(Yt,X)

    ARorange = [0.001 0.2:.2:.8];
    MAorange = [-.6:0.2:-0.2 0.01 .2:.2:.6];

    ARl = numel(ARorange); 
    MAl = numel(MAorange);

    LL = zeros(ARl,MAl);
    for ar_cnt = 1:ARl
        for ma_cnt = 1:MAl
            %disp(num2str([ARorange(ar_cnt) MAorange(ma_cnt)]))
            LL(ar_cnt,ma_cnt) = ARMA_SS_Lik(Yt,ARorange(ar_cnt),MAorange(ma_cnt),1,X);
        end
    end
    
    mLL_AFNI    = min(LL(:));
    [xs2,ys2]   = find(LL==mLL_AFNI);
    a           = ARorange(xs2);
    b           = MAorange(ys2);
    
    % what if there are more than one min?!
    a           = a(1); 
    b           = b(1);
    
end
