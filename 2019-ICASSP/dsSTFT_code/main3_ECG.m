close all ; clear all ;
% basic parameters for STFT
basicTF.win = 1201; %4096;
basicTF.hop = 21; %441;
basicTF.fs = 250;
basicTF.fr = 5/250; % frequency resolution/sampling freq
basicTF.feat = 'SST11'; % two option: STFT or SST11 (one window rejection)

% advanced parameters for STFT
advTF.num_tap = 1; % Number of tap in ConceFT
advTF.win_type = 'Gauss'; % Only 2-tap basis here (Hamming and its derivative)
advTF.Smo = 1; % alpha smoothing; 1 = no smoothing
advTF.Rej = 0; % The bandwidth of window rejection; 
advTF.ths = 1E-9; % Global threshold of STFT
advTF.HighFreq = 10/250; % highest frequency/sampling freq
advTF.LowFreq = 0.1/250; % lowest frequency/sampling freq
advTF.lpc = 0;
P.num_s = 1; % (harmonic 的根數, 只做簡單探討的話設1就好)
P.num_c = 1;

% parameters for cepstral representation
cepR.g = 0.3; % for generalized cepstrum: 設0 -> log cepstrum; 設2 -> autocorrelation (for single component 可以設大一點 (e.g., 2); for multiple component 經驗值: 0.1~0.3)
cepR.Tc=0; %1E-4; % Global threshold of cepstrum

x=importdata('ECG.mat');
x0 = x(:,2);
[a,b] = findpeaks(x0,'MinPeakHeight',0.4) ;
[irr] = interp1((b(2:end)+b(1:end-1))/2, 1000./(b(2:end)-b(1:end-1)), 1:length(x), 'cubic') ;

Trend = zeros(size(x0)) ;
for ii = 1: length(Trend)
	idx = [max(1,ii-50):min(length(Trend),ii+50)] ;
	Trend(ii) = median(x0(idx)) ;
end

x = x0 - Trend ;

x = resample(x,1,4);
t = [1:length(x)]'/basicTF.fs ;

fs = basicTF.fs;
% [b,a] = butter(9,0.01,'high');
% x = filter(b,a,x);
% y=resample(x(:,1),1,10);
[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(x, basicTF, advTF, cepR, P);
    



time_stamp = basicTF.hop/basicTF.fs;
figure(1)
subplot(231)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('STFT');

subplot(232)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
xlabel('time (s)'); ylabel('freqnency (s)');

subplot(233)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.999); axis xy; colormap(1-gray); axis([-inf inf 0 4]) ;
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Synchrosqueezed cepstrum');

subplot(234)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, abs(ceps), 0.999); axis xy; colormap(1-gray); axis([-inf inf 0 4]) ;
xlabel('time (s)'); ylabel('quefrency (sec)'); %title('Remapped cepstrum');

subplot(235)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tceps, 0.999); axis xy; colormap(1-gray); axis([-inf inf 0 4]) ;
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Remapped cepstrum');

subplot(236)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.999); axis xy; colormap(1-gray); axis([-inf inf 0 4]) ;
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Synchrosqueezed cepstrum');
hold on; plot(t, irr(1:4:end), 'r') ;
