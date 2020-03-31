function [Xp,Yp] = DrawMeSpectrum(Y,TR,dflag)

if ~exist('TR','var'); TR=1; end;
if ~exist('dflag','var'); dflag=0; end;

T = numel(Y);
Y = Y./std(Y);
Fs = 1./TR; 
n = 2^nextpow2(T);
f = Fs*(0:(n/2))/n;
FF = fft(Y,n);
P = abs(FF/n);

Xp=f;
Yp=P(1:n/2+1);

if dflag
    figure; hold on; box on; grid on; 
    plot(f,P(1:n/2+1)) 
end

end