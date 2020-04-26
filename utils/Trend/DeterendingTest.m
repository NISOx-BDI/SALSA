clear

SubID     = 'A00027167';
SesID     = 'DS2'; 
disp('=======================================')
PATH2AUX='/Users/sorooshafyouni/Home/GitClone/FILM2';
addpath([PATH2AUX '/utils/Trend'])
addpath('/Users/sorooshafyouni/Home/matlab/spm12')
%12487
disp('=====SET UP PATHS =============================')
%Raw Images (MMP feat output)
Path2ImgRaw=[PATH2AUX '/ExampleData/R.mpp'];
%Path2ImgDir = [Path2ImgRaw '/sub-' SubID '/ses-' SesID '/sub-' SubID '_ses-' SesID '_task-rest_acq-' num2str(TR*1000) '_bold.mpp'];
Path2ImgDir = ['/Users/sorooshafyouni/Home/GitClone/FILM2/Externals/R_test/sub-' SubID '_ses-' SesID '_task-rest_acq-645_bold_mpp'];
Path2Img    = [Path2ImgDir '/prefiltered_func_data_bet.nii'];







T = 500; 

t  = 1:T;
TR = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make 1/f noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% copied from: https://uk.mathworks.com/matlabcentral/answers/344786-generation-of-1-f-noise-using-matlab
fv = linspace(0, 1, 20);                                % Normalised Frequencies
a = 1./(1 + fv*2);                                      % Amplitudes Of ?1/f?
b = firls(42, fv, a);                                   % Filter Numerator Coefficients                                     % Filter Bode Plot
ns = rand(1, T);
invfn = filtfilt(b, 1, ns);                             % Create ?1/f? Noise

Y = invfn; %+ LowFreq(1:T);
Y = Y - mean(Y); % just demean so we can see the results easier.  
% similar to AFNI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pnum = 1+floor(TR.*T/150);
pnum = 3;
Yt2 = multpolyfit(repmat(1:T,1,1),Y,T,pnum);

%Just to check whether xt3 matches the matlab functions
Yt4 = Y - polyval(polyfit(t,Y,pnum),t);

% FSL what does %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filter = hp_fsl(T,100,TR);
Yft0   = filter*Y';


% SPM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SPMdefaultHPCutOff = 128;
% built-in SPM
%K = struct('RT',TR,'row',1:T,'HParam',SPMdefaultHPCutOff);
%Ydct = spm_filter(K,Y'); %DCT is happenening withing the function 
% Similar:
DCT_pnum  = fix(2*(T*TR)/SPMdefaultHPCutOff + 1);
DCT_hpf   = spm_dctmtx(T,DCT_pnum);
dctY = Y' - DCT_hpf*(DCT_hpf'*Y');


% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on; 
plot(t,Y,'LineWidth',1.5,'LineStyle','-.')
plot(t,Yt2,'LineWidth',1.5,'LineStyle','-.','Color','k')
plot(t,Yft0,'LineWidth',1.5,'LineStyle','-.','Color','b')
plot(t,Yt4)
plot(t,dctY,'LineWidth',1.5,'LineStyle','-.','Color','r')
legend('raw','poly','hpf','polyfitmatlabfunc','DCT')
%plot(t,x-y)

% MSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean((HighFreq'-dctY).^2)
% mean((HighFreq-Yft0').^2)
% mean((HighFreq-Yt2).^2)
