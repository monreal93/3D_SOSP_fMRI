function [f ,X, period] = frequencyAnalysis(x,dt)
% return Fourier transform amplitudes of gradient waveform x
% with raster time dt
% [dt] = s
n = length(x);
t = (0:1:n−1)*dt;
% convert waveforms to desired units
% x = mr.convert(x,'Hz/m','mT/m');
% calculate periods of waveform
[pks,locs] = findpeaks(x,t);
period = diff(locs);
fs = 1/dt;
% zero−centered frequency range
f = (−n/2:n/2−1)*(fs/n);
% FFT amplitude divided by sequence duration
X = abs(fftshift((fft(x))))/(t(end));
end
