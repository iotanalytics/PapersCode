clear all ; close all ;

% basic parameters for STFT
basicTF.win = 9001; %4096;
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
P.num_s = 1; % (harmonic 的根數, 只做簡單探討的話設1就好)
P.num_c = 1;

% parameters for cepstral representation
cepR.g = 0.3; % for generalized cepstrum: 設0 -> log cepstrum; 設2 -> autocorrelation (for single component 可以設大一點 (e.g., 2); for multiple component 經驗值: 0.1~0.3)
cepR.Tc=0; % Global threshold of cepstrum

load('TBME2013-PPGRR-Benchmark_R3/data/0333_8min.mat');
load('TBME2013-PPGRR-Benchmark_R3/data/0148_8min.mat');
load('TBME2013-PPGRR-Benchmark_R3/data/0104_8min.mat');
%load('TBME2013-PPGRR-Benchmark_R3/data/0370_8min.mat');
%load('TBME2013-PPGRR-Benchmark_R3/data/0009_8min.mat');
%load('TBME2013-PPGRR-Benchmark_R3/data/0016_8min.mat');
%load('TBME2013-PPGRR-Benchmark_R3/data/0329_8min.mat');
[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(signal.pleth.y - mean(signal.pleth.y), basicTF, advTF, cepR, P);

	% 121:400 means the frequency range 0.6Hz to 2.5 Hz.
	% 0.6 Hz means 36 heart beats per minute, which 
	% is kind of the lower limtation of cardiac activity
[cE] = CurveExt_M(abs(rtfr(121:500,:))', .2); 
cE = cE + 120 ;
IHR = tfrtic(cE)*basicTF.fs ;
time_stamp = basicTF.hop/basicTF.fs;
deshapeTime = 0:time_stamp:time_stamp*(size(rtfr,2)-1) ;

	% downsample the signal to 15 Hz
pleth_down = resample(signal.pleth.y - mean(signal.pleth.y), 1, 20) ;
[~, ~, tfrsq, ~, tfrsqtic] = ConceFT_STFT(pleth_down, 0, 0.5, 1e-3, 1, 15*30+1, 1, 6, 1, 0, 0, 0) ;
ConceftTime = [1:length(pleth_down)] / (basicTF.fs/20) ;


	% convert the IHR estimated from de-shape 
	% to that determined by SST, for the purpose of reconstruction
alpha0 = (tfrsqtic(2)-tfrsqtic(1)) * (basicTF.fs/20) ;
c = round(interp1(deshapeTime, IHR, [1:length(pleth_down)]/(basicTF.fs/20), 'linear', 'extrap') ./ alpha0) ;


HDpart = [] ;
RESPpart = [] ;
	% reconstruction from the frequency band ranging from 
	% 0.2Hz below the estimated IHR to the end.
RR = ceil(0.2/ (basicTF.fs/20)) ;
alpha = tfrsqtic(2)-tfrsqtic(1) ;
for kk = 1: length(c)
    idx = max(1,c(kk)-RR)+1: size(tfrsq, 1) ;
    HDpart(kk) = 2*sum(tfrsq(idx,kk),1)*((basicTF.fs/20)/2/0.5)*(alpha)/(basicTF.fs/20) ;
    idx = 1:max(1,c(kk)-RR) ;
    RESPpart(kk) = 2*sum(tfrsq(idx,kk),1)*((basicTF.fs/20)/2/0.5)*(alpha)/(basicTF.fs/20) ;
end

	% then run de-shape again on the RESP part, and you will get the correct frequency.
basicTF2 = basicTF ;
basicTF2.win = basicTF2.fs * 6 * 6 + 1 ;
[tfr2, ceps2, tceps2, tfrr2, rtfr2, tfrsq2, tfrtic] = CFPH(upsample(RESPpart,20,1)', basicTF2, advTF, cepR, P);




fs = basicTF.fs;
dt = basicTF.hop/basicTF.fs;



subplot(121)
plot(resample(signal.pleth.y - mean(signal.pleth.y),1,20))
hold on; plot(real(RESPpart),'r','linewidth',2)
axis([1000 1500 -inf inf]) ;

subplot(122)
plot(resample(signal.pleth.y - mean(signal.pleth.y),1,20))
hold on; plot(real(HDpart),'r','linewidth',2)
axis([1000 1500 -inf inf]) ;



[tfrC, cepsC, tcepsC, tfrrC, rtfrC, tfrsqC, tfrticC] = CFPH(signal.co2.y - mean(signal.co2.y), basicTF, advTF, cepR, P);


figure ; 
time_stamp = basicTF.hop/basicTF.fs;
subplot(131) ;
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrticC*basicTF.fs, tfrrC, 0.999); axis xy; colormap(1-gray); axis([-inf inf 0 1]) ;

subplot(132) ; 
imageSQ(deshapeTime, tfrtic*basicTF.fs, abs(tfrr), 0.98); axis([-inf inf 0 1])

subplot(133) ; 
imageSQ(deshapeTime, tfrtic*basicTF.fs, abs(tfrr2), 0.98); axis([-inf inf 0 1])
