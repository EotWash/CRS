function [amp_spect_den, freq] = asd2(time_series, Ts, smooth_width, poly_fit_terms, win_pointer, varargin)
% ASD2  calculate SMOOTHED, WINDOWED amplitude spectral density of a time series
% 
%    [amp_spect_den, freq] = asd2(time_series, Ts, smooth_width, poly_fit_terms, window_pointer, window_options)
%
%    time_series    the time series, 
%    Ts             the sampling time in seconds (time between consectutive samples)
%    smooth_width   number of raw fft bins to average across (e.g. 9) -
%                   this is like the number of averages. MUST BE ODD!
%    poly_fit_terms order of polynomial to be removed from the data,
%                   (we use polyfit) 0 = DC, 1 = DC and best fit line, 
%                   2 = DC, line, and parabola, etc.  (1 is good)
%    window_pointer    pointer to the built-in window, (e.g. @hann)
%                      see 'help window' for more information on windows 
%    window_options    options for the windows, per 'help window'
%  
%    amp_spect_den  the amplitude spectral density of the signal, and
%    freq           an optional output parameter, the corresonding freq vector
%
%    e.g. for data set called st2y_time_series, which was saved at 1/2048 time
%    steps, you might type:
%    [st2y_asd, freq] = asd2(st2y_time_series, 1/2048, 9, 1, @hann);
%
%    Note - the final bin width will be
%      BW = number of averages / length of time series in sec. 
%      e.g. a 1000 sec time series with 50 averages will return an
%      asd with bin width (freq resolution) of 50/(1000 sec) = 50 mHz
%
%    You can also call asd2 with only the first 2 or 3 or 4 inputs, e.g.     
%      [amp_spect_den, freq] = asd2(time_series, Ts)
%      [amp_spect_den, freq] = asd2(time_series, Ts, smooth_width)
%      [amp_spect_den, freq] = asd2(time_series, Ts, smooth_width, poly_fit_terms)
%
%    You can all call asd2 with [] as the inputs for:
%       smooth_width,  poly_fit_terms, and win_pointer, e.g.
%
%    [amp_spect_den, freq] = asd2(time_series, Ts, [], [], win_pointer)
%
%    inputs which are [] or not included will be set to the default values:
%    smooth_width = 9, poly_fit_terms = 1, and win_pointer = @hann.
%
%       e.g.
%  
%    [amp_spect_den, freq] = asd2(time_series, Ts, [], 3); 
%             is equivalent to 
%    [amp_spect_den, freq] = asd2(time_series, Ts,  9, 3, @hann);
%   
%    if time_series has units of xx, then amp_spect_den has units of xx-rms/rtHz
%    (assuming that Ts is in seconds)
%    
%    The first ASD is computed according to the formula
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
%    THEN, adjacent bins (smooth_width) are rms averaged together
%    BTL wrote asd.m on Oct 8, 2001
%
%    WINDOWS - the default window is the hanning window (@hann)
%    and the default option we use for the hanning window is 'periodic'
%    All windows are normalized so that the expected average power of a stationary 
%    signal will not change. i.e. NOT the Default matlab windows have a max amp of 1.
%    The scaling is done per:  
%     normalized_window = sqrt(1/mean(raw_window.^2)) * raw_window;
%       (so the mean of the normalized_window.^2 is 1)
%       
%   To use a DIFFERENT WINDOW, call with the window pointer to the window
%   you want (e.g. the tukey window is called with @tukeywin)
%   and add the options, as defined in >> help window, 
%     e.g.
%   [amp_spect_den, freq] = asd2(time_series, Ts,  9, 3, @tukeywin, 0.2);
%      or 
%   [amp_spect_den, freq] = asd2(time_series, Ts,  9, 3, @hann, 'symmetric');
%      see test_asd2.m for some other examples
%
%   WHAT'S THE BIG DEAL? - we average across freq bins, rather than
%   averaging sequential time series. This reduces the impact of the
%   shenanigans which result from windows, drifts and offsets. This can
%   be quite important if the spectra are not white. 
%
%   see also tfe2 and coh2 and mccs2
%
%    BTL wrote asd2.m on Sept 27, 2012
% the averaging does funny stuff to the DC term - don't trust it.
%
%  BTL update on March 5 to allow user options for windows
% 
%           IMPORTANT  2017  UPDATE
%
% Conor and BTL have adjusted the row/column direction of the outputs.
% If the input time series is a row, then freq and amp_spect_den are rows
% If the input is a column, then both outputs are columns.
% The old version had the freq vector as a row, requiring lots of freq.'
%     in analysis. this should now be just freq. 
%     this is a pain, because all of your old code will break.
%     To run old code, you can replace asd2 with asd2_pre20170929
%


% we do several steps:
% 0) Check inputs, set defaults for non-existant or empty values.
% 1) detrend the time series
%     by default we remove a 1st order poly, i.e. mean and slope,
%     but the user can pick something else.
% 2) apply a single window to the whole time series, 
%    periodic hann window by default (0 at one end, one point away from 0 at the
%    other), user can pick. The window is normalized so that the mean POWER
%    (amplitude squared) of the window is 1. 
% 3) take the ASD of the detrended, windowed, data
% 4) average the ASD by adjecent bins.


%% 0) setup stuff
debugging = false;

% make the input into a column vector

if ~exist('time_series','var') || isempty(time_series)
    error('input time_series not defined, see help asd2')
end

if ~exist('Ts','var') || isempty(Ts)
    error('input Ts not defined, see help asd2')
end

[rows, cols] = size(time_series);
if (rows > 1) && (cols > 1)
    error('input time series must be a vector')
end

if (cols > 1)    % input is a row vector, change it to a column
    input_was_a_row = true;
    time_series = time_series.';  % flip it
else
    input_was_a_row = false;
end

if debugging == true
    disp('debugging mode is on');
    
    if input_was_a_row
        disp('I flipped the input')
    else
        disp('I didn''t flip the input')
    end
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

% if ~exist('poly_fit_terms','var') || isempty(poly_fit_terms)
%     poly_fit_terms = 1;
% end

if ~exist('win_pointer','var') || isempty(win_pointer)
    win_pointer = @hann;
end

if isempty(varargin)  % could add checking and processing here.
    user_window_args = false;
else
    user_window_args = true;
end


%% 1) do the detrending
orig_pnts      = length(time_series);

if ~exist('poly_fit_terms','var') || isempty(poly_fit_terms)
    
    detrended_data = detrend(time_series); % Detrending via inbuilt function

else                                       % Or detrend with Polyfit, as prev.
    time           = Ts * (1:orig_pnts)';

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

    fit_coefs      = polyfit(time(1:step:end), time_series(1:step:end), poly_fit_terms);
    detrended_data = time_series - polyval(fit_coefs, time);
end

if debugging == true
    figure
    pp = plot(time, time_series, 'b', time, detrended_data, 'm');
    set(pp,'LineWidth',2)
    title('orig data vs. detrended data')
    xlabel('time (sec)')
    ylabel('mag')
    legend('orig data', 'detrended data')
    grid on
    FillPage('w')
    IDfig

    disp(' ')
    disp([' using a ',num2str(poly_fit_terms),' order fit'])
    disp('poly fit coefs are: ')
    disp(fit_coefs)
    
end

clear  time_series

%% 2) apply normalized window
if user_window_args == true   % make the window using user's options
    if debugging == true
        disp('using user window args')
        disp(varargin(:))
    end
    
    % make the window
    win = window(win_pointer, orig_pnts, varargin{:});
        
else                          % make the window using default or preset values
    
    if strcmp(func2str(win_pointer), 'hann')  % replaces strcmp(win_pointer, '@hann')
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
end   % end basic window constuction

% normalize the window so that the average signal.^2 (i.e. POWER) does not change.
win_norm = sqrt(1/mean(win.^2));
win = win * win_norm;   

detrended_windowed_data = win .* detrended_data;

if debugging == true
    figure
    subplot(211)
    pp = plot(time, detrended_data, 'b', time, detrended_windowed_data, 'm');
    set(pp,'LineWidth',2)
    title('detrended data vs. windowed, detrended data')
    xlabel('time (sec)')
    ylabel('mag')
    legend('detrended data','windowed and detrended data')
    grid on
    
    subplot(212)
    pp = plot(time, win);
    set(pp,'LineWidth',2)
    title('Window')
    xlabel('time (sec)')
    ylabel('mag')
    legend(['window type: ',func2str(win_pointer)])
    grid on
    
    FillPage('t')
    IDfig
    
    pow_detrend_series = mean(detrended_data.^2);
    pow_win_detrend_series = mean(detrended_windowed_data.^2);
    disp('avg power:')
    disp(['before win = ',num2str(pow_detrend_series)]);
    disp([' after win = ',num2str(pow_win_detrend_series)]);
    
end

clear time raw_win detrended_data;

%% 3) take asd
[full_asd, full_freq] = asd(detrended_windowed_data, Ts);

clear detrended_windowed_data
%% 4) average over adjacent bins
%avg_rule = 'single_DC_term';
avg_rule = 'uniform_freq_vector';


if strcmp(avg_rule, 'single_DC_term')
    % f(DC)            -> F(DC)
    % f(2)..f(1+width) -> F(2)
    %
    % so the number of new points is
    % 1 (for the DC term) + floor((old_fft_points - 1(for the DC term))/width)
    
    old_fft_points    = length(full_asd);
    new_fft_points    = 1 + floor((old_fft_points-1)/smooth_width);
    last_point_to_use = 1 + (new_fft_points-1) * smooth_width;
    
    amp_spect_den = zeros(new_fft_points,1);
    
    amp_spect_den(1) = full_asd(1);  % mirror the DC term (which should be zero)
    amp_spect_den(2:end) = sqrt(mean(reshape(...
        full_asd(2:last_point_to_use).^2, (new_fft_points-1), smooth_width), 2));

    % we make a new, rectangular matrix,
    % 1 column per final output value, and 'smooth_width' elements in each column,
    % then rms average down each column.
    
    freq = zeros(size(amp_spect_den));
    freq(1) = 0;
    freq(2:end) = sqrt(mean(reshape(...
        full_freq(2:last_point_to_use).^2, smooth_width, (new_fft_points-1)), 1));

end

if strcmp(avg_rule, 'uniform_freq_vector')
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
        
        old_fft_points    = length(full_asd);
        new_fft_points    = 1 + floor((old_fft_points-edge)/smooth_width);
        % note - the floor() ensures the averaging matrix is rectangular
        % and fully populated. We will probably not use a few high freq
        % points from the original full_fft.
        last_point_to_use = edge + (new_fft_points-1) * smooth_width;
        
        if (edge + smooth_width) > old_fft_points 
            error('your smooth_width is too large (or you data is too short)')
        end
        
        amp_spect_den = zeros(new_fft_points,1);
        
        if edge>1
            % assume the old DC term is 0, but average the non-zero terms into the
            % new DC term (odd, I admit)
            amp_spect_den(1) = sqrt(sum(full_asd(1:edge).^2)/(edge-1));
        else
            amp_spect_den(1) = full_asd(1);  % we are not averaging, and don't want to div by 0
        end
        
        amp_spect_den(2:end) = sqrt(mean(...
            reshape(full_asd((edge+1):last_point_to_use).^2, smooth_width, (new_fft_points-1)), 1));
    % we make a new, rectangular matrix,
    % 1 column per final output value, and 'smooth_width' elements in each column,
    % then rms average down each column.
        
        dF = full_freq(smooth_width+1);  % plus 1 because freq vector, DC is term 1
        freq = dF * (0:(new_fft_points-1)).'; % this is now a column vector
        
        % the asd is now a column vector.
        % if the input was a row vector, make the outputs row vectors also
        if input_was_a_row
            amp_spect_den = amp_spect_den.';
            freq = freq.';
        end
        warning('asd2:flipped',['\nThe asd2 freq output is now typically a column vector\n',...
            ' it used to be a row vector\n\n',...
            ' This may break old code, please be careful \n', ...
            ' see help asd2 for details \n',...
            '   this message only appears once \n'])
        warning('OFF', 'asd2:flipped')
    else
        error('right now, smooth_width must be odd')
    end
    
end


if debugging == true
    figure
    loglog(full_freq, full_asd, 'b', freq, amp_spect_den,'m')
    title('full ASD vs. averaged asd')
    xlabel('freq (Hz)')
    ylabel('mag / rtHz')
    legend('orig','avg''ed')
    
    warning('ON', 'asd2:flipped')

end




    dc_ave = floor((smooth_width-1)/2);  % minus 1 because DC term is element 1;

    old_fft_points    = length(full_asd);
    new_fft_points    = 1 + floor((old_fft_points-edge)/smooth_width);
    % note - the floor() ensures the averaging matrix is rectangular
    % and fully populated. We will probably not use a few high freq
    % points from the original full_fft.
    last_point_to_use = edge + (new_fft_points-1) * smooth_width;

    if (edge + smooth_width) > old_fft_points 
        error('your smooth_width is too large (or you data is too short)')
    end

    amp_spect_den = zeros(new_fft_points,1);

    if edge>1
        % assume the old DC term is 0, but average the non-zero terms into the
        % new DC term (odd, I admit)
        amp_spect_den(1) = sqrt(sum(full_asd(1:edge).^2)/(edge-1));
    else
        amp_spect_den(1) = full_asd(1);  % we are not averaging, and don't want to div by 0
    end

    amp_spect_den(2:end) = sqrt(mean(...
        reshape(full_asd((edge+1):last_point_to_use).^2, smooth_width, (new_fft_points-1)), 1));
% we make a new, rectangular matrix,
% 1 column per final output value, and 'smooth_width' elements in each column,
% then rms average down each column.

    dF = full_freq(smooth_width+1);  % plus 1 because freq vector, DC is term 1
    freq = dF * (0:(new_fft_points-1)).'; % this is now a column vector


