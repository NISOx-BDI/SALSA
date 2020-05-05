clear
% Parameters and signal
Fs = 1000; 
dt = 1/Fs; 
t = 0:dt:1-dt;
x = cos(2*pi*100*t)';
N = length(x);

% Add a bit of noise
x = x + 0.001*randn(N, 1);

[freq,psdx] = dft_periodogram(x,Fs);
[pxx, ~] = periodogram(x, rectwin(N), N, Fs);
[sdp,wp] = periodogram(x,[],N,Fs);

figure; 
hold on; 
plot(freq,log10(psdx))
plot(freq,log10(pxx))


function [freq,hpsdx] = dft_periodogram(x,Fs)

N = length(x);
% Use window function of your choice. Rectangular window means no window.
win = rectwin(N); 
% The actual number of points used when calculating the DFT might be higher
% this can help in speeding up the computation (next power of 2), but won't
% increase the frequency resolution (only frequency spacing is affected).
% See: https://dsp.stackexchange.com/a/32161/8202
NFFT = N; % N=2^20
xw   = x.*win; % Apply windowing here (periodogram does it on it's own)
xdft = fft(xw, NFFT);

psdx = abs(xdft).^2; %xdft.*conj(xdft); % psdx = abs(xdft).^2
hpsdx = psdx(1:NFFT/2+1); % Take only first half of the PSD
% Calculate the scaling factor. 
% Note that in case of rectangular window this is simply N.
S    = sum(win.^2);
SF   = (2/(Fs*S));
hpsdx = SF * hpsdx;
freq = 0:Fs/NFFT:Fs/2;
% Calculate the corresponding PSD using MATLAB's in-built
[pxx, ~] = periodogram(x, rectwin(N), NFFT, Fs);

phase0 = phase(xdft);
% size(phase)
% ts0 = real(ifft( sqrt(hpsdx/SF) .* exp(1j*phase0),NFFT));
% ts0(1:5)
% Plot the PSD's
% plot(freq, 10*log10(psdx), 'r')
% grid on
% hold on
% plot(freq, 10*log10(pxx), 'b--')
% title('PSD calculation comparison')
% xlabel('Frequency [Hz]')
% ylabel('Power/Frequency [dB/Hz]')
% legend({'DFT', 'periodogram'})

end