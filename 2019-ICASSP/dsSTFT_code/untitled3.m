%close all ;clear all ;

%clear
%close all



%load sig

sig=ppg_syn.data(1).ppg.v;

%%

% basic parameters for STFT
basicTF.win = 1001; %4096;
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
advTF.HighFreq = 10/basicTF.fs; % highest frequency/sampling freq
advTF.LowFreq = 0.1/basicTF.fs; % lowest frequency/sampling freq
advTF.lpc = 0;
P.num_s = 2; % (harmonic ���ڼ�, �u��²�汴�Q���ܳ]1�N�n)
P.num_c = 1;
% parameters for cepstral representation
cepR.g = 0.5; % for generalized cepstrum: �]0 -> log cepstrum; �]2 -> autocorrelation (for single component �i�H�]�j�@�I (e.g., 2); for multiple component �g���: 0.1~0.3)
cepR.Tc = 0; % Global threshold of cepstrum

%x = importdata('PhyCam005_Walk_Walk_100Hz_pulse.mat');
%x = importdata( 'test_sig.mat');
x = [sig',sig'];

x = x(1:2*floor(end/2),:) ;
z0 = x(:,1) - mean(x(:,1)) ;
%[x, trend] = deTrend([1:length(z0)]'/basicTF.fs, z0, 0.2) ;
x = detrend(z0) ;

[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(x, basicTF, advTF, cepR, P);
fs = basicTF.fs;
dt = basicTF.hop/basicTF.fs;

time_stamp = basicTF.hop/basicTF.fs;


%%
%[cR] = CurveExt_M(abs(rtfr(1:200,:))', 2);

[cR2] = CurveExt_M(abs(rtfr(201:400,:))', 2)+200;

[cR3] = CurveExt_M(abs(rtfr(401:600,:))', 2)+400;

%%
 load PhyCam005_Walk_Walk_ECG_IBI_edit.txt
 irr = interp1(PhyCam005_Walk_Walk_ECG_IBI_edit(:,1), 1000./ PhyCam005_Walk_Walk_ECG_IBI_edit(:,2), 0:time_stamp:time_stamp*(size(rtfr,2)-1), 'cubic','extrap');


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
plot(0:1/basicTF.fs:1/basicTF.fs*(size(x)-1),x,'k','linewidth',1.2); axis([0 inf -12 12]) %colormap(1-gray);
set(gca,'xtick',[],'ytick',[])
%xlabel('time (s)'); ylabel('amplitude'); title('PPG');

subplot(222)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.995); axis xy; colormap(1-gray);axis([0 inf 0 4])
set(gca,'xtick',[],'ytick',[])
%xlabel('time (s)'); ylabel('frequency (Hz)'); title('STFT');

% subplot(222)
% %imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
% imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrr, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 4]) ;
% xlabel('time (s)'); ylabel('freqnency (Hz)'); title('Remapped cepstrum');%title('SST');

subplot(223)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 4]);
set(gca,'xtick',[],'ytick',[])
%xlabel('time (s)'); ylabel('frequency (Hz)');  title('SST');%title('Remapped cepstrum');

subplot(224)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 4]) ;
set(gca,'xtick',[],'ytick',[])
%xlabel('time (s)'); ylabel('frequency (Hz)'); title('de-shape SST');%title('Synchrosqueezed cepstrum');
%hold on; %plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), irr,'r')


%%
figure

subplot(221)
plot(0:1/basicTF.fs:1/basicTF.fs*(size(x)-1),x,'k','linewidth',1.2); %axis([0 inf -12 12]) %colormap(1-gray);
set(gca,'xtick',[],'ytick',[])
%xlabel('time (s)'); ylabel('amplitude'); title('PPG');

subplot(222)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.995); axis xy; colormap(1-gray);%axis([0 inf 0 6])
set(gca,'xtick',[],'ytick',[])
%xlabel('time (s)'); ylabel('frequency (Hz)'); title('STFT');

% subplot(222)
% %imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
% imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrr, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 4]) ;
% xlabel('time (s)'); ylabel('freqnency (Hz)'); title('Remapped cepstrum');%title('SST');

subplot(223)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy; colormap(1-gray); %axis([0 inf 0 6]);
set(gca,'xtick',[],'ytick',[])
%xlabel('time (s)'); ylabel('frequency (Hz)');  title('SST');%title('Remapped cepstrum');
%plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), irr,'r')


subplot(224)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.995); axis xy; colormap(1-gray); %axis([0 inf 0 6]) ;
set(gca,'xtick',[],'ytick',[])
%xlabel('time (s)'); ylabel('frequency (Hz)'); title('de-shape SST');%title('Synchrosqueezed cepstrum');
%hold on; plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic(cR)*fs,'r')
hold on; plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic(cR3)*fs,'r')


%%

%%
figure

subplot(221)
plot(0:1/basicTF.fs:1/basicTF.fs*(size(x)-1),x,'k','linewidth',1.2); axis([0 inf -12 12]) %colormap(1-gray);
set(gca,'xtick',[],'ytick',[])
%xlabel('time (s)'); ylabel('amplitude'); title('PPG');

subplot(222)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.995); axis xy; colormap(1-gray);axis([0 inf 0 6])
set(gca,'xtick',[],'ytick',[])
%xlabel('time (s)'); ylabel('frequency (Hz)'); title('STFT');

% subplot(222)
% %imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
% imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrr, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 4]) ;
% xlabel('time (s)'); ylabel('freqnency (Hz)'); title('Remapped cepstrum');%title('SST');

subplot(223)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 6]);
set(gca,'xtick',[],'ytick',[])
%xlabel('time (s)'); ylabel('frequency (Hz)');  title('SST');%title('Remapped cepstrum');
%plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), irr,'r')


subplot(224)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.995); axis xy; colormap(1-gray); axis([0 inf 0 6]) ;
set(gca,'xtick',[],'ytick',[])
%xlabel('time (s)'); ylabel('frequency (Hz)'); title('de-shape SST');%title('Synchrosqueezed cepstrum');
hold on; plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic(cR)*fs,'r')
plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic(cR2)*fs,'k')
plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic(cR3)*fs,'b')

%%

mean(tfrtic(cR)*fs)
std(tfrtic(cR)*fs)

%%

mean(tfrtic(cR2)*fs)
std(tfrtic(cR2)*fs)

mean(tfrtic(cR3)*fs)
std(tfrtic(cR3)*fs)

%%

load('0009_8min.mat');
%load('TBME2013-PPGRR-Benchmark_R3/data/0329_8min.mat');
[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(signal.pleth.y - mean(signal.pleth.y), basicTF, advTF, cepR, P);

[tfrC, cepsC, tcepsC, tfrrC, rtfrC, tfrsqC, tfrticC] = CFPH(signal.co2.y - mean(signal.co2.y), basicTF, advTF, cepR, P);

Trend = zeros(size(signal.ecg.y)) ;
for ii = 1: length(Trend)
    idx = [max(1,ii-15):min(length(Trend),ii+15)] ;
    Trend(ii) = median(signal.ecg.y(idx)) ;
end

x = signal.ecg.y - Trend ;

[tfrE, cepsE, tcepsE, tfrrE, rtfrE, tfrsqE, tfrticE] = CFPH(x, basicTF, advTF, cepR, P);
fs = basicTF.fs;
dt = basicTF.hop/basicTF.fs;

time_stamp = basicTF.hop/basicTF.fs;

%%


[cR] = CurveExt_M(abs(rtfr(1:200,:))', 2);
[cE] = CurveExt_M(abs(rtfr(201:400,:))', .2); cE = cE + 200 ;

% advTF.Rej = 2;
% advTF.num_tap = 10;
%%

figure(1)
subplot(231)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');

subplot(232)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
xlabel('time (s)'); ylabel('freqnency (s)');

subplot(234)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.999); axis xy; colormap(1-gray); axis([-inf inf 0 4]) ;
xlabel('time (s)'); ylabel('frequency (Hz)');

subplot(235)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.999); axis xy; colormap(1-gray); axis([-inf inf 0 2]) 
hold on ; plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic(cR)*fs*5,'b');
plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic(cR)*fs,'r');
plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic(cE)*fs,'r');
xlabel('time (s)'); ylabel('frequency (Hz)');

subplot(233)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfrC, 0.999); axis xy; colormap(1-gray); axis([-inf inf 0 2]) ;
hold on; plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic(cR)*fs,'r');
xlabel('time (s)'); ylabel('frequency (Hz)'); 

subplot(236)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfrE, 0.999); axis xy; colormap(1-gray); axis([-inf inf 0 2])
hold on; plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic(cE)*fs,'r');
xlabel('time (s)'); ylabel('frequency (Hz)');
