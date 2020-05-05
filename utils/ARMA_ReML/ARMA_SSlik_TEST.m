clear

T = 500; 
addpath('/Users/sorooshafyouni/Downloads/dlm-master')
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
aNom  = [1 -0.8];%ARcoeffs
b  = [1 0.4];%MAcoeffs
e  = randn(T,nRlz); % generate gaussian white noise
Ytall = filter(b,aNom,e); % generate y

%Ytall = randn(T,100); 

Ytall = Ytall./std(Ytall);
Ytall = Ytall - mean(Ytall); 

% set a 2D grid -----------------------------------
ARorange = [0.001 0.2:.2:.8];
MAorange = [-.6:0.2:-0.2 0.01 .2:.2:.6];

for iRlz = 1    
   Yt = Ytall(:,iRlz);
   %--------------------
   Yt = Yt - mean(Yt); 
   Yt = Yt./std(Yt);
    
    ARl = numel(ARorange); 
    MAl = numel(MAorange);   

    LL_AFNI = zeros(ARl,MAl);
    for ar_cnt = 1:ARl
        for ma_cnt = 1:MAl    
            l(ar_cnt,ma_cnt) = armalik(Yt,ARorange(ar_cnt),MAorange(ma_cnt),1);
        end
    end
    
   ml = min(l(:));
   [xl,yl] = find(l==ml);
   pp(:,iRlz) = [ARorange(xl) MAorange(yl)]
end