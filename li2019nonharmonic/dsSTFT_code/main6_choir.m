close all ; clear all ;


%addpath('D:\Audio_word_toolbox\toolbox\MIRtoolbox1.3.3\MIRToolbox');
% basic parameters for STFT
basicTF.win = 2205; %4096;
basicTF.hop = 441; %441;
basicTF.fs = 44100;
basicTF.fr = 5;
basicTF.feat = 'SST11'; % STFT or SST11
% advanced parameters for STFT
advTF.num_tap = 1;%num_tap(cc);
advTF.win_type = 'Gauss'; %{'Gaussian','Thomson','multipeak','SWCE'%};
advTF.Smo = 1;
advTF.Rej = 0;
advTF.ths = 1E-9;
advTF.HighFreq = 1000/44100;
advTF.LowFreq = 50/44100;
advTF.lpc = 0;
P.num_s = 1;
P.num_c = 1;
% parameters for cepstral representation
cepR.g = 0.3; % g(cc);%0.06;
cepR.Tc = 0;


% [x fs]=wavread('clarinet1.wav');
% [x fs]=wavread('C:\Users\user\Dropbox\中研院合作音樂家003\DATA\Violin sonata\VS05_Mozart\VS05_Mozart_audio.wav');
[x fs]=audioread('CH01_Bach_audio.wav');
% [x fs]=mp3read('C:\Users\user\Dropbox\中研院合作音樂家003\midi+mp3 cut\lieder\5 brahms\mp3 cut.mp3');
% x = mean(x,2); x = x(1:44100*6);% x = x(44100*6.5:44100*11.55);
% wavwrite(x,'CH01_Bach_audio.wav');
[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(x, basicTF, advTF, cepR, P);

fs = basicTF.fs;
dt = basicTF.hop/basicTF.fs;
time_stamp = basicTF.hop/basicTF.fs;


figure(1)
subplot(221)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr.^0.5, 0.9999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); %title('STFT');

subplot(222)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq.^0.5, 0.9999); axis xy;
xlabel('time (s)'); ylabel('freqnency (s)');

subplot(223)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr.^0.5, 0.9999); axis xy; colormap(1-gray); axis([0 inf 0 800]) ;
xlabel('time (s)'); ylabel('frequency (Hz)'); 

subplot(224)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr.^0.5, 0.9999); axis xy; colormap(1-gray); axis([0 time_stamp*(size(rtfr,2)-1) 0 800]) ;
xlabel('time (s)'); ylabel('frequency (Hz)'); 

load('CH01_Bach_tutti.mat');
hold on ;
    for ci = 1:length(GTNotes)
        z=round(440.*2.^((GTNotes{ci}(2,:)-69)./12));%./basicTF.fr)+1;
        plot(time_stamp*GTNotes{ci}(1,1):time_stamp:time_stamp*GTNotes{ci}(1,end), z,'r','linewidth', 2);
    end
hold off;





scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/4]) ;
plot([1:length(x)]/fs, x,'k')
set(gca,'fontsize', 18) ; xlabel('Time (s)') ; ylabel('arbitrary unit') ;
axis([4.2 4.7 -inf inf]) ;
