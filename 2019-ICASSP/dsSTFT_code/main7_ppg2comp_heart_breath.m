clear all ; close all ;

% basic parameters for STFT
basicTF.win = 4001; %4096;
basicTF.hop = 301; %441;
basicTF.fs = 300;
basicTF.fr = 0.005; % frequency resolution/sampling freq
basicTF.feat = 'SST11'; % two option: STFT or SST11 (one window rejection)

% advanced parameters for STFT
advTF.num_tap = 1; % Number of tap in ConceFT
advTF.win_type = 'Gauss'; % Only 2-tap basis here (Hamming and its derivative)
advTF.Smo = 1; % alpha smoothing; 1 = no smoothing
advTF.Rej = 0; % The bandwidth of window rejection; 
advTF.ths = 1E-9; % Global threshold of STFT
advTF.HighFreq = 10/300; % highest frequency/sampling freq
advTF.LowFreq = 0.05/300; % lowest frequency/sampling freq
advTF.lpc = 0;
P.num_s = 1; % (harmonic ���ڼ�, �u��²�汴�Q���ܳ]1�N�n)
P.num_c = 1;

% parameters for cepstral representation
cepR.g = 0.3; % for generalized cepstrum: �]0 -> log cepstrum; �]2 -> autocorrelation (for single component �i�H�]�j�@�I (e.g., 2); for multiple component �g���: 0.1~0.3)
cepR.Tc=0; % Global threshold of cepstrum

%load('TBME2013-PPGRR-Benchmark_R3/data/0333_8min.mat');
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


[cR] = CurveExt_M(abs(rtfrE(601:1000,:))', 1);cR = cR+600;
[cE] = CurveExt_M(abs(rtfrE(301:600,:))', .2); cE = cE + 300 ;

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

%%

mean(tfrtic(cE)*fs)
std(tfrtic(cE)*fs)

% mean(tfrtic(cR)*fs*5)
% std(tfrtic(cR)*fs*5)


%%

figure
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfrE, 0.999); axis xy; colormap(1-gray); axis([-inf inf 0 6])
%hold on; plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrticE(cE)*fs,'r');
%plot(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrticE(cR)*fs,'r');
%xlabel('time (s)'); ylabel('frequency (Hz)');
set(gca,'xtick',[],'ytick',[])


mean(tfrtic(cE)*fs)
std(tfrtic(cE)*fs)

mean( tfrticE(cR)*fs)
std( tfrticE(cR)*fs)
