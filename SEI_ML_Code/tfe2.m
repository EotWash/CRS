function [ tf_estimate, freq ] = tfe2( time_series1, time_series2 , Ts, smooth_width, poly_fit_terms, win_pointer)
%tfe  calculate SMOOTHED, WINDOWED transfer function between 2 time series
% 
%    [tf_estimate, freq] = tfe2(time_series1, time_series2, Ts, smooth_width, poly_fit_terms, win_pointer)
%
%     - required inputs-
%    time_series1   the 2 time series, must be the same length
%    time_series2
%    Ts             the sampling time in seconds (time between consectutive samples)
%     - optional inputs -
%    smooth_width   number of raw fft bins to average across (default is 9) -
%                   this is like the number of averages. MUST BE ODD!
%    poly_fit_terms order of polynomial to be removed from the data,
%                   (we use polyfit) 0 = DC, 1 = DC and best fit line, 
%                   2 = DC, line, and parabola, etc.  (default is 1)
%    win_pointer    pointer to the built-in window, (default is @hann)
%                     see 'help window' for more information on windows 
%  
%     - outputs -
%    tf_estimate   estimate of the transfer function: series 2/ series 1
%    freq           an optional output parameter, the corresonding freq vector
%       Note - the final bin width will be
%         BW = number of averages / length of time series in sec. 
%         e.g. a 1000 sec time series with 50 averages will return a
%         tfe with bin width (freq resolution) of 50/(1000 sec) = 50 mHz
%
%    e.g. for data sets called st1z_time_series and st2z_time_series, 
%         which were saved at 1/2048 time steps, you might type:
%    [stg2_o_stg1_TF, freq] = tfe2(st1z_time_series, st2z_time_series, 1/2048, 9, 1, @hann);
%
%    You can also call tfe2 with only the first 2 or 3 or 4 inputs, e.g.     
%      [tf_estimate, freq] = tfe2(time_series1, time_series2, Ts)
%      [tf_estimate, freq] = tfe2(time_series1, time_series2, Ts, smooth_width)
%      [tf_estimate, freq] = tfe2(time_series1, time_series2, Ts, smooth_width, poly_fit_terms)
%
%    You can all call tfe2 with [] as the inputs for:
%       smooth_width,  poly_fit_terms, and win_pointer, e.g.
%
%    [tf_estimate, freq] = tfe2(time_series1, time_series2, Ts, [], [], win_pointer)
%
%    inputs which are [] or not included will be set to the default values:
%    smooth_width = 9, poly_fit_terms = 1, and win_pointer = @hann.
%       e.g.
%    [tf_estimate, freq] = tfe2(time_series1, time_series2, Ts, [], 3); 
%       is equivalent to 
%    [tf_estimate, freq] = tfe2(time_series1, time_series2, Ts,  9, 3, @hann);
%   
%    tf_estimate has dimesions of (signal 2 dimension)/ (signal 1 dimension) 
%       it is defined as TF12 = (<P21>) ./ <P22>)
%       where Pxy = (1/N)^2 * (1/BW) * 2 * fft(time_seriesX) .* conj(fft(time_seriesY)) 
%       ie the single sided power spectrum or cross spectrum, as appropriate.
%       <Pxy> is the result of averaging (mean) adjacent bins of the
%       spectrum. The averaging is done before the division.
%    
%    N is the number of points in the time series, and 
%    BW is the bandwidth of a freq bin (ie the bin width), BW = 1/(N*Ts);
%
%    Except for the first and last points,
%    Matlabs ffts are double sided, and the PSDs are single sided, hence the 2.
%
%    The DC term and the Nyquist freq term have the 2, which is bogus, and the averaging
%    there is weird. Don't believe these points. (Besides, the bias is
%    removed by the detrending, but will be partly re-introduced by the
%    window. Then you average in the first (smooth_ramp-1)/2 spectral
%    points resulting in ??. Treat this term with extreme case. It shows up
%    at freq = 0 so doesn't plot for loglog or semilogx plots.
%
%   WHAT'S THE BIG DEAL? - we average across freq bins, rather than
%   averaging sequential time series. This reduces the impact of the
%   shenanigans which result from windows, drifts and offsets. This can
%   be quite important if the spectra are not white. 
%
%    BTL wrote asd2.m on Sept 27, 2012
%    and adapted it to coh2.m on Oct 16, 2012
%    and adapted it to tfe2.m on Oct 29, 2012

% steps are:
% 0) Check inputs, set defaults for non-existant or empty values.
% 1) detrend each time series
%     by default we remove a 1st order poly, i.e. mean and slope,
%     but the user can pick something else.
% 2) apply a single window to each time series, 
%    periodic hann window by default (0 at one end, one point away from 0 at the
%    other), user can pick. The window is normalized so that the mean POWER
%    (amplitude squared) of the window is 1. 
% 3) calculate PXX and PYX of each detrended, windowed, data
% 4) average PXX, PYX by adjecent bins into pxx and pyx
%    use averaged spec to calc cxy.

%% 0) check inputs

debugging = false;

if ~exist('time_series1','var') || isempty(time_series1)
    error('input time_series1 not defined, see help tfe2')
end

if ~exist('time_series2','var') || isempty(time_series2)
    error('input time_series2 not defined, see help tfe2')
end


% make the input into a column vector
[rows, cols1] = size(time_series1);
if (rows > 1) && (cols1 > 1)
    error('input time series 1 must be a vector')
end

[rows, cols2] = size(time_series2);
if (rows > 1) && (cols2 > 1)
    error('input time series 2 must be a vector')
end

if (cols1 > 1)    % input1 is a row vector, change it
    time_series1 = time_series1.';  % flip it
end
if (cols2 > 1)    % input2 is a row vector, change it
    time_series2 = time_series2.';  % flip it
end

if length(time_series1) ~= length(time_series2)
    error('input time series must be the same length')
end


if ~exist('Ts','var') || isempty(Ts)
    error('input Ts not defined, see help tfe2')
end


if ~exist('smooth_width','var') || isempty(smooth_width)
    smooth_width = 9;
end


if smooth_width ~= round(smooth_width)
    disp('smooth width must be an integer')
    smooth_width = ceil(smooth_width);
    disp(['RESETTING smooth_width to ',num2str(smooth_width)])
end

if smooth_width <1
    error('smooth width must be >= 1')
end

if ~exist('poly_fit_terms','var') || isempty(poly_fit_terms)
    poly_fit_terms = 1;
end

if ~exist('win_pointer','var') || isempty(win_pointer)
    win_pointer = @hann;
end


%% 1) detrend the data
orig_pnts       = length(time_series1);
time            = Ts * (1:orig_pnts)';
 

 if orig_pnts < 1e4
    step = 1;
  elseif orig_pnts < 1e5
    step = 10;
  elseif orig_pnts < 1e6
    step = 100;
  elseif orig_pnts < 1e7
    step = 1E3;
  elseif  orig_pnts < 1e8
    step = 1E4;
   else 
    step = 1E5;
    cprintf([1 0 0.5],'You are insane, how much data do you think that we can handle?\n')
 end
    

fit_coefs1      = polyfit(time(1:step:end), time_series1(1:step:end), poly_fit_terms);
detrended_data1 = time_series1 - polyval(fit_coefs1, time);

fit_coefs2      = polyfit(time(1:step:end), time_series2(1:step:end), poly_fit_terms);
detrended_data2 = time_series2 - polyval(fit_coefs2, time);

if debugging == true
    figure
    subplot(211)
    pp = plot(time, time_series1, 'b', time, detrended_data1, 'm');
    set(pp,'LineWidth',2)
    title('orig data 1 vs. detrended data 1')
    xlabel('time (sec)')
    ylabel('mag')
    legend('orig data', 'detrended data')
    grid on
    
    subplot(212)
    pp = plot(time, time_series2, 'b', time, detrended_data2, 'm');
    set(pp,'LineWidth',2)
    title('orig data 2 vs. detrended data 2')
    xlabel('time (sec)')
    ylabel('mag')
    legend('orig data', 'detrended data')
    grid on
    FillPage('t')
    IDfig

    disp(' ')
    disp([' using a ',num2str(poly_fit_terms),' order fit'])
    disp('poly fit coefs for time series 1 are: ')
    disp(fit_coefs1)
    disp('poly fit coefs for time series 2 are: ')
    disp(fit_coefs2)
    
end

clear  time_series1 time_series2

%% 2) apply normalized window

if strcmp(win_pointer, '@hann')
    win_args = true;   % are there any arguments?
    win_opts = 'periodic';
else
    win_args = false;
    win_opts = [];
end

if win_args == true
    win = window(win_pointer, orig_pnts, win_opts);
else
    win = window(win_pointer, orig_pnts);
end

win_norm = sqrt(1/mean(win.^2));
win = win * win_norm;

detrended_windowed_data1 = win .* detrended_data1;
detrended_windowed_data2 = win .* detrended_data2;

if debugging == true
    figure
    subplot(211)
    pp = plot(time, detrended_data1, 'b', time, detrended_windowed_data1, 'm');
    set(pp,'LineWidth',2)
    title('detrended data 1 vs. windowed, detrended data 1')
    xlabel('time (sec)')
    ylabel('mag')
    legend('detrended data','windowed and detrended data')
    grid on
    
    subplot(212)
    pp = plot(time, detrended_data2, 'b', time, detrended_windowed_data2, 'm');
    set(pp,'LineWidth',2)
    title('detrended data 2 vs. windowed, detrended data 2')
    xlabel('time (sec)')
    ylabel('mag')
    legend('detrended data','windowed and detrended data')
    grid on
    FillPage('t')
    IDfig

    pow_detrend_series1 = mean(detrended_data1.^2);
    pow_win_detrend_series1 = mean(detrended_windowed_data1.^2);
    pow_detrend_series2 = mean(detrended_data2.^2);
    pow_win_detrend_series2 = mean(detrended_windowed_data2.^2);

    disp('avg power:')
    disp(['series 1, before win = ',num2str(pow_detrend_series1)]);
    disp(['series 1, after win = ',num2str(pow_win_detrend_series1)]);
    disp(' ')
    disp(['series 2, before win = ',num2str(pow_detrend_series2)]);
    disp(['series 2, after win = ',num2str(pow_win_detrend_series2)]);
end

clear time win detrended_data1 detrended_data2;

%% 3) calc the psds and csds
fX_big = fft(detrended_windowed_data1);
fX     = fX_big(1:ceil(end/2));
clear fX_big

fY_big = fft(detrended_windowed_data2);
fY     = fY_big(1:ceil(end/2));
clear fY_big

N  = orig_pnts;
BW = 1/(N*Ts);

full_freq = (0:length(fX)-1)*BW;

% from BTLs asd.m code:
%N  = length(time_series);
%BW = 1/(N*Ts);
%sig_fft        = fft(time_series);
%amp_spect_temp = sqrt(2) * (1/N) * sqrt(sig_fft.* conj(sig_fft)) * (1/sqrt(BW));

% for the coherence and transfer functions,
% we DO NOT NEED TO do all the normalizations, because they all cancel
% 
%       HOWEVER
% we will DO THEM ANYWAYS, just to avoid copy-paste troubles in the future
% PXX is the single-sided, mean squared spectral density


PXX = 2 * (1/N)^2 * (1/BW) * (fX .* conj(fX));
PYX = 2 * (1/N)^2 * (1/BW) * (fY .* conj(fX));
clear fX fY



%% 3)averaging - use the 'constant frequency vector method

if smooth_width == 1 + 2*round((smooth_width-1)/2)  % is smooth_width odd?
    % let W = smooth_width, which is odd.
    % the DC term will be W wide, for a 2 sided FFT,
    % centered about DC term, ie
    % let edge = (W-1)/2, the DC goes from -edge -> + edge (for 2 sided ASD)
    % or DC -> +edge (for our single sided ASD)
    % f(DC:edge)               -> F(DC)
    % f(edge+1)..f(edge+width) -> F(2)
    %
    % so the number of new points is
    % 1 (for the DC term) + floor((old_fft_points - 1(for the DC term))/width)
    edge = (smooth_width-1)/2 + 1;  % plus 1 because DC term is element 1;
    
    old_fft_points    = length(PXX);
    new_fft_points    = 1 + floor((old_fft_points-edge)/smooth_width);
    % note - the floor() ensures the averaging matrix is rectangular
    % and fully populated. We will probably not use a few high freq
    % points from the original full_fft.
    last_point_to_use = edge + (new_fft_points-1) * smooth_width;
    
    if (edge + smooth_width) > old_fft_points
        error('your smooth_width is too large (or you data is too short)')
    end
    
    % these are the averaged versions
    tf_estimate = zeros(new_fft_points,1);
    pxx = zeros(size(tf_estimate));
    pyx = zeros(size(tf_estimate));
    
    % average the old DC term and nearby bits to new DC term
    % this term is crap.
    pxx(1) = mean(PXX(1:edge));
    pyx(1) = mean(PYX(1:edge));
    
    pxx(2:end) = mean(...
        reshape(PXX((edge+1):last_point_to_use), smooth_width, (new_fft_points-1)), 1);
    clear PXX
    pyx(2:end) = mean(...
        reshape(PYX((edge+1):last_point_to_use), smooth_width, (new_fft_points-1)), 1);
    clear PYX
    % we make a new, rectangular matrix,
    % 1 column per final output value, and 'smooth_width' elements in each column,
    % then average down each column.
    tf_estimate = pyx ./ pxx ;
    
    dF = full_freq(smooth_width+1);  % plus 1 because freq vector, DC is term 1
    freq = dF * (0:(new_fft_points-1));
else
    error('right now, smooth_width must be odd')
end



end

