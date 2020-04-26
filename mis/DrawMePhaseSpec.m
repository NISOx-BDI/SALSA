function [Xp,Yp,theta] = DrawMePhaseSpec(Y,TR,dflag)
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
f   = Fs*(0:(T/2))/T;
FF  = fft(Y)/T;
P   = abs(FF).^2; % normalise by the number of paddings // P = abs(FF);

Xp  = f;
Yp  = P(1:T/2+1,:);

theta = unwrap(angle(FF));
theta = theta(1:T/2+1,:);

%--- spectrum flattness
% nn       = T/2+1;
% gmean    = exp(sum(log(Yp(2:end,:)))./nn);
% specflat = gmean./mean(Yp(2:end,:));

if dflag
    figure; hold on; box on; grid on; 
    plot(f,theta) 
end

end