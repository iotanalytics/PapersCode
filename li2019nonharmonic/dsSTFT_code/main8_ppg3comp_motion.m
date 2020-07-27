% basic parameters for STFT
basicTF.win = 2001; %4096;
basicTF.hop = 101; %441;
basicTF.fs = 125;
basicTF.fr = 0.02; % frequency resolution
basicTF.feat = 'SST11'; % two option: STFT or SST11 (one window rejection)

% advanced parameters for STFT
advTF.num_tap = 1; % Number of tap in ConceFT
advTF.win_type = 'Gauss'; % Only 2-tap basis here (Hamming and its derivative)
advTF.Smo = 1; % alpha smoothing; 1 = no smoothing
advTF.Rej = 0; % The bandwidth of window rejection; 
advTF.ths = 1E-9; % Global threshold of STFT
advTF.lpc = 0;
advTF.HighFreq = 10/125; % highest frequency/sampling freq
advTF.LowFreq = 0.05/125; % lowest frequency/sampling freq
P.num_s = 1; % (harmonic 的根數, 只做簡單探討的話設1就好)
P.num_c = 1;
% parameters for cepstral representation
cepR.g = 0.3; % for generalized cepstrum: 設0 -> log cepstrum; 設2 -> autocorrelation (for single component 可以設大一點 (e.g., 2); for multiple component 經驗值: 0.1~0.3)
cepR.Tc= 0; % Global threshold of cepstrum

load('DATA_02_TYPE02.mat');
[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(sig(2,:)- mean(sig(2,:)), basicTF, advTF, cepR, P);
fs = basicTF.fs;
dt = basicTF.hop/basicTF.fs;

time_stamp = basicTF.hop/basicTF.fs;


figure(1)
subplot(221)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.995); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('STFT');

subplot(222)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
xlabel('time (s)'); ylabel('freqnency (s)'); %title('SST');

subplot(223)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tceps, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 6])
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Remapped cepstrum');

subplot(224)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 6])
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Synchrosqueezed cepstrum');
