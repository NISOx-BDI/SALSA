function [Xp,Yp,specflat] = DrawMeSpectrum(Y,TR,dflag)
% Gets the spectrum of given signal Y
% Y should be TxV
% TR: repeatition time [default: 1]
% dfalg: if switched, it plots the spectrum
%
%
% SA, Ox, 2020

if ~exist('TR','var'); TR=1; end;
if ~exist('dflag','var'); dflag=0; end;

T   = size(Y,1);
Y   = Y./std(Y);
Fs  = 1./TR; 
%orig ===
% n   = 2^nextpow2(T); % this will be done in fft
% f   = Fs*(0:(n/2))/n;
% FF  = fft(Y,n);
%==
Xp   = Fs*(0:round(T/2))/T; % that is freq range
FF  = fft(Y);
P   = abs(FF).^2/T; % normalise by the number of paddings // P = abs(FF);

Yp  = P(1:round(T/2)+1,:);

%--- spectrum flattness
nn       = T/2+1;
gmean    = exp(sum(log(Yp(2:end,:)))./nn);
specflat = gmean./mean(Yp(2:end,:));

if dflag==1
    figure; hold on; box on; grid on; 
    plot(Xp,P(1:round(T/2)+1)) 
elseif dflag==-1
    figure; 
    hold on; box on; grid on;
    magntd = abs(FF);
    plot(Xp,magntd(1:T/2+1))     
end

end