% basic parameters for STFT
basicTF.win = 3001; %4096;
basicTF.hop = 101; %441;
basicTF.fs = 100;
basicTF.fr = 0.5/100; % frequency resolution/sampling freq
basicTF.feat = 'SST11'; % two option: STFT or SST11 (one window rejection)

% advanced parameters for STFT
advTF.num_tap = 1; % Number of tap in ConceFT
advTF.win_type = 'Gauss'; % Only 2-tap basis here (Hamming and its derivative)
advTF.Smo = 1; % alpha smoothing; 1 = no smoothing
advTF.Rej = 0; % The bandwidth of window rejection; 
advTF.ths = 1E-9; % Global threshold of STFT
advTF.HighFreq = 2/100; % highest frequency/sampling freq
advTF.LowFreq = 0.02/100; % lowest frequency/sampling freq
advTF.lpc = 0;
P.num_s = 1; % (harmonic 的根數, 只做簡單探討的話設1就好)
P.num_c = 1;

% parameters for cepstral representation
cepR.g = 0.3; % for generalized cepstrum: 設0 -> log cepstrum; 設2 -> autocorrelation (for single component 建議設大一點 (e.g., 2); for multiple component 經驗值: 0.1~0.3)
cepR.Tc= 0; % Global threshold of cepstrum

load('Su.mat');
[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(x1, basicTF, advTF, cepR, P);
t = [1:length(x1)]'/basicTF.fs ;
t = t(1:basicTF.hop:end) ;

% advTF.Rej = 2;
% advTF.num_tap = 10;
% [~, ~, ~, rtfr2, ~] = CFPH(x1, basicTF, advTF, cepR, P);
% basicTF.feat = 'SST11';
% [~, ~, ~, rtfr3, ~] = CFPH(x1, basicTF, advTF, cepR, P);
fs = basicTF.fs;
dt = basicTF.hop/basicTF.fs;

figure(1)
subplot(221)
imageSQ(t, tfrtic*basicTF.fs, tfr, 0.995); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('STFT');

subplot(222)
imageSQ(t, tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
xlabel('time (s)'); ylabel('freqnency (s)'); %title('SST');

subplot(223)
imageSQ(t, tfrtic*basicTF.fs, tceps, 0.995); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Remapped cepstrum');

subplot(224)
imageSQ(t, tfrtic*basicTF.fs, rtfr, 0.995); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Synchrosqueezed cepstrum');
