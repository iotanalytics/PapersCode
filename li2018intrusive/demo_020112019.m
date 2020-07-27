

clear 
close all
clc

samplingrate = 100;

fftlength = 10000;

sig = csvread('periodicalRepiration1.csv');

sig = sig(2000:11000);

sig_bp = bandpass(sig-mean(sig),[0.7 8],samplingrate);

%%

sig_bp3 = bandpass(sig-mean(sig),[0.6 2],samplingrate);

[Pxx, F] = pyulear(sig_bp3,4,10000,samplingrate);
figure
plot(F,Pxx)
xlabel('Frequency (Hz)')
title('Autoregressive (AR) Spectral Analysis')
xlim([0 5])

%%

sig_bp = bandpass(sig-mean(sig),[0.7 8],samplingrate);

     [yu,yl] = envelope(sig_bp,60,'peak');
     
     %%band pass the respiratory only
     

     

    [rrpks,rrlocs] = findpeaks(yl,'MinPeakWidth',50);
    
    
    
    figure
subplot(211)
plot(0:1/samplingrate:(length(sig)-1)/samplingrate,sig_bp);
xlabel('Time (s)')
ylabel('Amplitude')
title 'Vibration Signal'
xlim([0 30])

 subplot(212)
 hold on
plot(0:1/samplingrate:(length(sig)-1)/samplingrate,sig_bp);
plot(0:1/samplingrate:(length(sig)-1)/samplingrate,yu);
plot(0:1/samplingrate:(length(sig)-1)/samplingrate,yl);
%plot(0:1/samplingrate:(length(sig)-1)/samplingrate,imf(:,1));
hold off
xlabel('Time (s)')
ylabel('Amplitude')
title 'Envelopes'
xlim([0 30])
box on
    
syl = fft(yl-mean(yl),fftlength);
Psyl = syl.*conj(syl)/fftlength;

syu = fft(yu-mean(yu),fftlength);
Psyu = syu.*conj(syu)/fftlength;



%%          

%windd = hann(64);
%yu_hann = 

      %    yu_bp=bandpass(yl,[0.1,0.7],samplingrate);
     windd = 300;
      yu = yu(windd:end-windd);
      yl = yu(windd:end-windd);
     
yu_hp=bandpass(yu,[0.1 0.6],samplingrate);
yl_hp=bandpass(yl,[0.1 0.6],samplingrate);

order =6;
          
    % [imf, residual, info] = emd(yu_bp,'Display',0);
[Pxx1, F] = pyulear(yu_hp,order,10000,samplingrate);
[Pxx2, F] = pyulear(yl_hp,order,10000,samplingrate);
[Pxx3, F] = pyulear([yl_hp;yu_hp],order,10000,samplingrate);

%Pxx3 = Pxx1/max(Pxx1)+Pxx2/max(Pxx2);

figure
hold on
plot(F,Pxx1/max(Pxx1))
plot(F,Pxx2/max(Pxx2))
plot(F,Pxx3/max(Pxx3))
xlabel('Frequency (Hz)')
title('Autoregressive (AR) Spectral Analysis')
legend('upper','lower','aver')
xlim([0 1])
     