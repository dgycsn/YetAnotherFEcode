% MYFUNCTION
%
% [freq, spectrummag, spectrumangle] = fft_yFs(y, Fs)
%
% computes one-sided spectrum given a history y and a sampling frequency Fs
% as inputs.

function [freq, spectrummag, spectrumangle] = fft_yFs(y, Fs)
N = length(y);
spectrum = fft(y);
spectrummag(1)=abs(spectrum(1))/N;
spectrummag(2:N/2+1)=abs(spectrum(2:N/2+1))/N*2;
spectrumangle(1:N/2+1)=angle(spectrum(1:N/2+1));
freq = Fs/N*(0:N/2);
end