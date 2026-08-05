function result = analyze_frequency(trace_matrix,profile,spec)
%ANALYZE_FREQUENCY Compute direct FFT and representative CWT summaries.
% Role and frame rate are read from the complete channel profile.
% Inputs:
%   trace_matrix Numeric [frames x ROI] matrix.
%   profile Scalar bound profile requiring role and positive frame_rate Hz.
%   spec Scalar struct requiring min_freq_hz, max_freq_hz,
%        wavelet_voices_per_octave, and wavelet_name.
% Output:
%   result Struct containing per-ROI FFT summaries and one representative
%          CWT. Nonfinite samples are linearly filled before analysis.
trace_matrix = double(trace_matrix);
if ~ismatrix(trace_matrix) || isempty(trace_matrix)
    error('AAA:Functions:InvalidFrequencyTrace','Trace matrix must be nonempty.');
end
frame_rate = double(profile.frame_rate);
max_hz = min(double(spec.max_freq_hz),frame_rate/2-eps);
min_hz = double(spec.min_freq_hz);
nrois = size(trace_matrix,2);
frequency=[]; amplitude=[]; peak_hz=NaN(nrois,1); peak_amp=NaN(nrois,1);
signal_rms=NaN(nrois,1);
for idx=1:nrois
    x=sanitize(trace_matrix(:,idx)); signal_rms(idx)=rms(x);
    [f,a]=fft_spectrum(x,frame_rate);
    if isempty(frequency), frequency=f; amplitude=NaN(numel(f),nrois); end
    amplitude(:,idx)=a;
    valid=f>=min_hz & f<=max_hz;
    if any(valid)
        indices=find(valid); [peak_amp(idx),local]=max(a(valid));
        peak_hz(idx)=f(indices(local));
    end
end
[~,representative]=max(peak_amp,[],'omitnan');
if isempty(representative) || ~isfinite(representative), representative=1; end
x=sanitize(trace_matrix(:,representative));
wavelet=representative_cwt(x,frame_rate,spec);
result=struct('role',string(profile.role),'frame_rate',frame_rate, ...
    'parameters',spec,'time',(0:size(trace_matrix,1)-1)'/frame_rate, ...
    'trace_matrix',trace_matrix,'fft',struct('frequency',frequency,'amplitude',amplitude), ...
    'roi',struct('peak_frequency_hz',peak_hz,'peak_amplitude',peak_amp, ...
        'signal_rms',signal_rms), ...
    'representative',struct('roi_index',representative,'trace',x, ...
        'wavelet',wavelet), ...
    'summary',struct('nrois',nrois, ...
        'median_peak_frequency_hz',median(peak_hz,'omitnan'), ...
        'median_peak_amplitude',median(peak_amp,'omitnan')));
end
function x=sanitize(x)
x=double(x(:));
if all(~isfinite(x)), x=zeros(size(x)); else, x=fillmissing(x,'linear','EndValues','nearest'); end
x=x-mean(x,'omitnan');
end
function [f,a]=fft_spectrum(x,fs)
n=numel(x);
if n<8, f=zeros(0,1); a=zeros(0,1); return; end
y=fft(x); p=abs(y/n); a=p(1:floor(n/2)+1); if numel(a)>2, a(2:end-1)=2*a(2:end-1); end
f=fs*(0:floor(n/2))'/n;
end
function value=representative_cwt(x,fs,spec)
value=struct('available',false,'coefficients',[],'frequency',[], ...
    'time',(0:numel(x)-1)'/fs,'method',string(spec.wavelet_name));
if exist('cwt','file')~=2, return; end
try
    [coefficients,frequency]=cwt(x,char(string(spec.wavelet_name)),fs, ...
        'VoicesPerOctave',double(spec.wavelet_voices_per_octave));
    value.available=true; value.coefficients=coefficients; value.frequency=frequency;
catch ME
    value.error=struct('identifier',string(ME.identifier),'message',string(ME.message));
end
end
