clear
close all
clc

A = importdata('905H_log(1).txt');


gamma = A(:,1);
impedance = A(:,2);
reflectivity = A(:,3);
density = A(:,4);
pwave = A(:,5);
synthetic = A(:,6);
trace = A(:,7);

dt = 0.001;
dl = length(A);

figure
subplot(611)
plot(0:dt:dt*(dl-1),gamma);
xlim([0 dt*(dl-1)])
title 'gamma'
subplot(612)
plot(0:dt:dt*(dl-1),impedance);
title 'impedance'
xlim([0 dt*(dl-1)])
subplot(613)
plot(0:dt:dt*(dl-1),reflectivity)
title 'reflectivity'
xlim([0 dt*(dl-1)])
subplot(614)
plot(0:dt:dt*(dl-1),density);
title 'density'
xlim([0 dt*(dl-1)])
subplot(615)
plot(0:dt:dt*(dl-1),pwave)
title 'P-wave'
xlim([0 dt*(dl-1)])
subplot(616)
plot(0:dt:dt*(dl-1),trace);
title 'seismic'
xlim([0 dt*(dl-1)])


%%

 sig = synthetic;
 sig2 = impedance;
 sig3 = reflectivity;
 sig4 = density;



% basic parameters for STFT
basicTF.win = 501; %4096;
basicTF.hop = 11; %441;
basicTF.fs = 500;
basicTF.fr = 0.01; % frequency resolution/sampling freq
basicTF.feat = 'SST11'; % two option: STFT or SST11 (one window rejection)

% advanced parameters for STFT
advTF.num_tap =1; % Number of tap in ConceFT
advTF.win_type = 'Gauss'; % Only 2-tap basis here (Hamming and its derivative)
advTF.Smo = 0.2; % alpha smoothing; 1 = no smoothing
advTF.Rej = 0; % The bandwidth of window rejection; 
advTF.ths = 1E-9; % Global threshold of STFT
advTF.HighFreq = 50/basicTF.fs; % highest frequency/sampling freq
advTF.LowFreq = 10/basicTF.fs; % lowest frequency/sampling freq
advTF.lpc = 0;
P.num_s = 2; % (harmonic ���ڼ�, �u��²�汴�Q���ܳ]1�N�n)
P.num_c = 1;
% parameters for cepstral representation
cepR.g = 0.2; % for generalized cepstrum: �]0 -> log cepstrum; �]2 -> autocorrelation (for single component �i�H�]�j�@�I (e.g., 2); for multiple component �g���: 0.1~0.3)
cepR.Tc = 0; % Global threshold of cepstrum

%x = importdata('PhyCam005_Walk_Walk_100Hz_pulse.mat');
%x = importdata( 'test_sig.mat');
x = [sig,sig];

x = x(1:2*floor(end/2),:) ;
z0 = x(:,1) - mean(x(:,1)) ;
%[x, trend] = deTrend([1:length(z0)]'/basicTF.fs, z0, 0.2) ;
x = z0 ;

[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(x, basicTF, advTF, cepR, P);
fs = basicTF.fs;
dt = basicTF.hop/basicTF.fs;

time_stamp = basicTF.hop/basicTF.fs;


% load PhyCam005_Walk_Walk_ECG_IBI_edit.txt
% irr = interp1(PhyCam005_Walk_Walk_ECG_IBI_edit(:,1), 1000./ PhyCam005_Walk_Walk_ECG_IBI_edit(:,2), 0:time_stamp:time_stamp*(size(rtfr,2)-1), 'cubic','extrap');


%%
% figure(1)
% subplot(221)
% imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.995); axis xy; colormap(1-gray);
% xlabel('time (s)'); ylabel('frequency (Hz)'); %title('STFT');
% 
% subplot(223)
% %imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
% imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 4]) ;
% xlabel('time (s)'); ylabel('freqnency (s)'); %title('SST');
% 
% subplot(222)
% imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy; colormap(1-gray); 
% xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Remapped cepstrum');
% 
% subplot(224)
% imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 4]) ;
% xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Synchrosqueezed cepstrum');
% hold on; plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), irr,'r')


%%
figure

subplot(221)
plot(0:1/basicTF.fs:1/basicTF.fs*(size(x)-1),x,'k','linewidth',1.5); %axis([0 inf -12 12]) %colormap(1-gray);
xlabel('time (s)'); ylabel('amplitude'); title('PPG');

subplot(222)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.995); axis xy; colormap(1-gray);%axis([0 inf 0 4])
xlabel('time (s)'); ylabel('frequency (Hz)'); title('STFT');

% subplot(222)
% %imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
% imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrr, 0.995); axis xy; colormap(1-gray); %axis([0 inf 0 4]) ;
% xlabel('time (s)'); ylabel('freqnency (Hz)'); title('Remapped cepstrum');%title('SST');

subplot(223)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy; colormap(1-gray);% axis([0 inf 0 4]);
xlabel('time (s)'); ylabel('frequency (Hz)');  title('SST');%title('Remapped cepstrum');

subplot(224)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.995); axis xy; colormap(1-gray); %axis([0 inf 0 4]) ;
xlabel('time (s)'); ylabel('frequency (Hz)'); title('de-shape SST');%title('Synchrosqueezed cepstrum');
hold on; %plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), irr,'r')


