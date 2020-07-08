close all ;clear all ;

% basic parameters for STFT
basicTF.win = 1001; %4096;
basicTF.hop = 21; %441;
basicTF.fs = 100;
basicTF.fr = 0.02; % frequency resolution/sampling freq
basicTF.feat = 'SST11'; % two option: STFT or SST11 (one window rejection)

% advanced parameters for STFT
advTF.num_tap = 1; % Number of tap in ConceFT
advTF.win_type = 'Gauss'; % Only 2-tap basis here (Hamming and its derivative)
advTF.Smo = 1; % alpha smoothing; 1 = no smoothing
advTF.Rej = 0; % The bandwidth of window rejection; 
advTF.ths = 1E-9; % Global threshold of STFT
advTF.HighFreq = 8/100; % highest frequency/sampling freq
advTF.LowFreq = 0.2/100; % lowest frequency/sampling freq
advTF.lpc = 0;
P.num_s = 1; % (harmonic ���ڼ�, �u��²�汴�Q���ܳ]1�N�n)
P.num_c = 1;
% parameters for cepstral representation
cepR.g = 0.3; % for generalized cepstrum: �]0 -> log cepstrum; �]2 -> autocorrelation (for single component �i�H�]�j�@�I (e.g., 2); for multiple component �g���: 0.1~0.3)
cepR.Tc = 0; % Global threshold of cepstrum

%x = importdata('PhyCam005_Walk_Walk_100Hz_pulse.mat');
x = importdata('sig_test.mat');
x = x(1:2*floor(end/2),:) ;
z0 = x(:,1) - mean(x(:,1)) ;
%[x, trend] = deTrend([1:length(z0)]'/basicTF.fs, z0, 0.2) ;
x = z0 ;

[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(x, basicTF, advTF, cepR, P);
fs = basicTF.fs;
dt = basicTF.hop/basicTF.fs;

time_stamp = basicTF.hop/basicTF.fs;


load PhyCam005_Walk_Walk_ECG_IBI_edit.txt
irr = interp1(PhyCam005_Walk_Walk_ECG_IBI_edit(:,1), 1000./ PhyCam005_Walk_Walk_ECG_IBI_edit(:,2), 0:time_stamp:time_stamp*(size(rtfr,2)-1), 'cubic','extrap');


%%
figure(1)
subplot(221)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.995); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('STFT');

subplot(223)
%imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 4]) ;
xlabel('time (s)'); ylabel('freqnency (s)'); %title('SST');

subplot(222)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy; colormap(1-gray); 
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Remapped cepstrum');

subplot(224)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 4]) ;
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Synchrosqueezed cepstrum');
hold on; plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), irr,'r')


%%
figure(2)

subplot(221)
plot(0:basicTF.fr/2:basicTF.fr/2*(size(x)-1),x,'k','linewidth',1.5); axis([0 inf -0.01 0.04]) %colormap(1-gray);
xlabel('time (s)'); ylabel('amplitude'); title('PPG');

subplot(222)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.995); axis xy; colormap(1-gray);axis([0 inf 0 4])
xlabel('time (s)'); ylabel('frequency (Hz)'); title('STFT');

% subplot(222)
% %imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
% imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrr, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 4]) ;
% xlabel('time (s)'); ylabel('freqnency (Hz)'); title('Remapped cepstrum');%title('SST');

subplot(223)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 4]);
xlabel('time (s)'); ylabel('frequency (Hz)');  title('SST');%title('Remapped cepstrum');

subplot(224)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 4]) ;
xlabel('time (s)'); ylabel('frequency (Hz)'); title('de-shape SST');%title('Synchrosqueezed cepstrum');
hold on; plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), irr,'r')

