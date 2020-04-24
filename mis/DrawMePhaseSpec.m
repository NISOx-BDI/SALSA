function [Xp,Yp,theta,specflat] = DrawMePhaseSpec(Y,TR,dflag)
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
n   = 2^nextpow2(T);
f   = Fs*(0:(n/2))/n;
FF  = fft(Y,n);
P   = abs(FF).^2/n; % normalise by the number of paddings // P = abs(FF);

Xp  = f;
Yp  = P(1:n/2+1,:);

theta = angle(FF);
theta = theta(1:n/2+1,:);

%--- spectrum flattness
nn       = n/2+1;
gmean    = exp(sum(log(Yp(2:end,:)))./nn);
specflat = gmean./mean(Yp(2:end,:));

if dflag
    figure; hold on; box on; grid on; 
    plot(f,theta) 
end

end