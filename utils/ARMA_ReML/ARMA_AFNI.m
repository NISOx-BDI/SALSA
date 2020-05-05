clear

T = 500; 

% Simulate Data ---------
% sMdl = arima(1,0,1);
% sMdl.Constant = 0;
% sMdl.Variance = 1;
% sMdl.AR = {0.6};
% sMdl.MA = {-0.2};
% Yt = simulate(sMdl,T);

% random noise
% Yt = randn(500,1); 
% Yt = Yt - mean(Yt); 
% Yt = Yt./std(Y);

nRlz = 1000; 

%simulate data using filter
aNom  = [1 -0.2];%ARcoeffs
b  = [1 0.6];%MAcoeffs
e  = randn(T,nRlz); % generate gaussian white noise
Ytall = filter(b,aNom,e); % generate y

%Ytall = randn(T,100); 

Ytall = Ytall./std(Ytall);
Ytall = Ytall - mean(Ytall); 

%[xx,yy] = DrawMeSpectrum(Ytall,0.645);

% figure; hold on; 
% plot(xx,mean(yy,2))
% plot(f,(abs(h)./T/2).^2)

for iRlz = 1:nRlz
    iRlz

    tic;
    
    Yt = Ytall(:,iRlz);
    %--------------------
    Yt = Yt - mean(Yt); 
    Yt = Yt./std(Yt);

% set a 2D grid -----------------------------------
    ARorange = [0:.2:.8];
    MAorange = [-.8:.2:.8];

    ARl = numel(ARorange); 
    MAl = numel(MAorange);

    LL_AFNI = zeros(ARl,MAl);
    for ar_cnt = 1:ARl
        for ma_cnt = 1:MAl
            %disp(num2str([ARorange(ar_cnt) MAorange(ma_cnt)]))
            LL_AFNI(ar_cnt,ma_cnt) = ARMA_ReML_AFNI(ARorange(ar_cnt),MAorange(ma_cnt),Yt,ones(T,1));
        end
    end
    allLL_AFNI(:,:,iRlz) = LL_AFNI;
    
    mLL_AFNI = min(LL_AFNI(:));
    [xs2_AFNI,ys2_AFNI]   = find(LL_AFNI==mLL_AFNI);
    AFNIp = [ARorange(xs2_AFNI) MAorange(ys2_AFNI)];
    AFNIp(:,iRlz)

    toc
end
