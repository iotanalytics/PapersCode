% addpath('D:\Audio_word_toolbox\toolbox\MIRtoolbox1.3.3\MIRToolbox');
% basic parameters for STFT
basicTF.win = 6175; %4096;
basicTF.hop = 441; %441;
basicTF.fs = 44100;
basicTF.fr = 4;
basicTF.feat = 'SST11'; % STFT or SST11
% advanced parameters for STFT
advTF.num_tap = 1;%num_tap(cc);
advTF.win_type = 'Blackman'; %{'Gaussian','Thomson','multipeak','SWCE'%};
advTF.Smo = 1;
advTF.Rej = 0;
advTF.ths = 1E-9;
advTF.HighFreq = 3000/44100;
advTF.LowFreq = 50/44100;
advTF.lpc = 0;
P.num_s = 1;
P.num_c = 1;
% parameters for cepstral representation
cepR.g = 0.3; % g(cc);%0.06;
cepR.Tc=0;


[x fs]=wavread('026.wav');
x = mean(x,2); x = x(1:44100*15);
[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(x, basicTF, advTF, cepR, P);
fs = basicTF.fs;
dt = basicTF.hop/basicTF.fs;
time_stamp = basicTF.hop/basicTF.fs;

% advTF.Rej = 2;
% advTF.num_tap = 10;


figure(1)
subplot(221)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.995); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('STFT');

subplot(222)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
xlabel('time (s)'); ylabel('freqnency (s)'); %title('SST');

subplot(223)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tceps, 0.995); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Remapped cepstrum');

subplot(224)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.995); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('Synchrosqueezed cepstrum');
