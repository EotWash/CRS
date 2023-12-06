function [amp_spect,freq]=asd(time_series,Ts)
% ASD  calculates the amplitude spectral density of a time series
%    [amp_spect_den,freq]=asd(time_series,Ts), where
%    time_series    the time series, 
%    Ts             the sampling time in seconds (time between consectutive samples)
%    amp_spect_den  the amplitude spectral density of the signal, and
%    freq           an optional output parameter, the corresonding freq vector
%
%    WARNING: there is no windowing, you must do it yourself.
%    see asd2.m for a more user-friendly function! (BTL Sept 2012)
% 
%    if time_series has units of xx, then amp_spect_den has units of xx-rms/rtHz
%    (assuming that Ts is in seconds)
%    
%    it is computed according to the formula
%    asd = sqrt(2) * 1/N * sqrt(sig_fft .* conj(sig_fft)) * 1/sqrt(BW)
%    where 
%    sig_fft = fft(time_series),
%    N is the number of points in the time series, and 
%    BW is the bandwidth of a freq bin (ie the bin width), BW = 1/(N*Ts);
%    except for the first and last points,
%    The ffts are double sided, and the ASD is single sided, hence the root 2.
%    sig_fft is length N, and ASD in length N/2+1.
%    The DC term and the Nyquist freq term do not have the sqrt(2)
%    since they do not have a conjugate term in the two sided FFT
%
%    BTL Oct 8, 2001
%
%  PS - Remember, when averaging together several ASDs, you must use
%       the rms average, not the linear average. Use the linear average
%       for power spectra.



N  = length(time_series);
BW = 1/(N*Ts);

sig_fft        = fft(time_series);
amp_spect_temp = sqrt(2) * (1/N) * sqrt(sig_fft.* conj(sig_fft)) * (1/sqrt(BW));

if N/2 == round(N/2)       % are there an even number of points in the time series?
	N2 = N/2 + 1;
	amp_spect         = zeros(1,N2);
	amp_spect(1)      = amp_spect_temp(1)/sqrt(2);
	amp_spect(2:N2-1) = amp_spect_temp(2:N2-1);
	amp_spect(N2)     = amp_spect_temp(N2)/sqrt(2);
else
	N2 = (N+1)/2;
	amp_spect         = zeros(1,N2);
	amp_spect(1)      = amp_spect_temp(1)/sqrt(2);
	amp_spect(2:N2)   = amp_spect_temp(2:N2);
end

if nargout == 2
    freq = (0:N2-1)*BW;
    disp(freq)
end
