function [asdOut, fOut] = asdSmooth(timeSeriesIn, ts, logStep, nAve)
% Make smoother ASDs
% 
% Spectal densities are increasing jaggy at high-frequencies when taken 
% with a standard linear frequency axis and fixed number of averages. 
% It's possible to stitch together sections with different resolution and 
% averaging, or to use the LTPDA Toolbox's log-spaced spectral density, 
% but this function should do something in between, retaining the 
% ASD2-style frequency-domain averaging (for superior spectral
% lobe suppression)  and effectively auto-stitching sections with 
% different smooth-width.
%
% [asdOut fOut] = asdSmooth(timeSeriesIn/asdIn, ts/fIn,unit, logStep, nAve)
%
% 'timeSeriesIn' is the input time series
% 
% 'ts' is the time-step of the data
% 
% Note that non-uniform frequency vectors make RMS calculations less
% accurate!
%
% 'logStep' (optional) is the mutiple of frequency in which to increase the 
% number of averages. Default is 6. At f = 6*min(fIn) the 'smooth_width'
% increases by a factor of nAve. At 100*min(fIn) it is a factor of nAve^2.
%
% 'nAve' (optional) the multiplicative factor to increase the smoothing by at
% each frequency step. Currently must be a positive, odd integer. Default is 3.
% 
% If you have an existing spectrum and want to smooth it, that's possible 
% by adding two more input arguments that essentially override the normal 
% asd2 process. The number of averages (and step) might need adjustment. 
% NOTE: the first two inputs are ignored in this case, recommended values 
% for logStep and nAve are 8 and 2 respectively.
% 
% 'asdIn' is the input ASD
% 'fIn' is the input frequency vector. Must be linear and start from zero.
%
% The spectrum will be truncated to avoid partially filling the final bin.
% 
% Conor Mow-Lowry 21/11/2017

%% Input variable checking and parsing

if ~exist('timeSeriesIn','var') || isempty(timeSeriesIn)
   error('Missing input timeseries/ASD. Seek help.')
end

if ~exist('ts','var') || isempty(ts)
    error('Missing time-step/frequency. Seek help.')
end

[rows, cols] = size(timeSeriesIn);
if (rows > 1) && (cols > 1)
    error('input time series must be an array')
end

% Step 1: Check if the user is sending frequency or time-domain data
if length(ts)>1
    % We're now using frequency domain data, so we'll set a flag to
    % indicate this so we don't take an asd of an asd....

    warning('asdSmooth:frequency',['\n Both inputs are vectors and asdSmooth will now \n',...
        ' treat the inputs as an ASD and frequency vector \n',...
        ' and perform RMS averaging over neighbouring bins. \n \n', ...
        ' If you intended to take the spectrum of a time-series \n',...
        ' please ensure the second input has a single element. \n',...
        ' This message only appears once. \n'])
    warning('OFF', 'asdSmooth:frequency')

    fIn=ts;
    asdIn=timeSeriesIn;
    
    % Default, 9-average asd2 spectra want less averaging than time-series.
    if ~exist('logStep','var') || isempty(logStep)
        logStep = 8;
    end
    if ~exist('nAve','var') || isempty(nAve)
        nAve = 2;
    end
    
    if length(fIn) ~= length(asdIn)
        error('Frequency and ASD vectors must be the same size')
    end
    
elseif length(ts)==1
    % Taking the ASD of the input time series with no averaging
    [asdIn, fIn] = asd2(timeSeriesIn, ts, 1);
    
    % Faster increase in # of averages for time-series inputs.
    if ~exist('logStep','var') || isempty(logStep)
        logStep = 9;
    end
    if ~exist('nAve','var') || isempty(nAve)
        nAve = 3;
    end
end


df = fIn(2) - fIn(1);
if any(round(100*diff(fIn)/df) ~= 100)
    error('asdSmooth only works with linear frequency vectors.')
elseif all([fIn(1) (fIn(1) - df)] ~= 0)
    error('Frequency vector should start from zero.')
end

if logStep<=nAve
   nAve = logStep -1;
   warning(['nAve should be less than logStep. \n',...
       'Setting nAve = logStep - 1'])
end



%% Breaking spectrum into chunks

% We want each chunk to fill it's frequency multiple with an even multiple
% of the number of averages.

y = floor(logStep/nAve)*nAve;


% To solve for the number of sections, we need (courtesy of Wolfram Alpha)
m = length(fIn)*(y-1)/y + 1;

% The number of chunks is then
z = log(m)/log(y); 

nStitch = ceil(z);      % The minimum number of discrete sections
n= 1:nStitch;           % A vector of those
nPts = y.^n(1:nStitch-1);   % Number of points in all the filled chunks
nFin = length(fIn) - sum(nPts);      % Initial number of points in final chunk
nSec = floor(nFin/nAve^(nStitch-1)); % Number of smooth-widths in final chunk
nPts(end+1) = nSec*nAve^(nStitch-1); % New number of points in final chunk
nSec = nPts./nAve.^(n-1);            % Number of smooth-sections in each chunk

fIn = fIn(1:sum(nPts));              % Truncating data
asdIn = asdIn(1:sum(nPts));

% Starting points of the chunks
b1 = [1 1 + cumsum(nPts(1:end-1))];

% Endpoints of chunks
b2 = cumsum(nPts);

%% Building outputs
asdOut = [];
fOut = [];
for k = 1:nStitch
    % Applying Brian's ASD2 'reshape' for averaging and RMS averaging
    fChunk = reshape(fIn(b1(k):b2(k)), nAve^(k-1), nSec(k), 1);
    asdChunk = reshape(asdIn(b1(k):b2(k)).^2, nAve^(k-1), nSec(k), 1);

    % Accumulating output ASD
    fOut = [fOut mean(fChunk,1)];
    asdOut = [asdOut sqrt(mean(asdChunk,1))];
end

if exist('rows','var')
    if rows >1      % Flipping the ASD to match its input
        fOut = fOut.';
        asdOut = asdOut.';
    end
end






