%% new HR and RR estimation
clear
close all
clc

fs= 100;

%% read data
sig = csvread('joseHB.csv');
sigj = sig(4000:10000);

%% FIR filter to remove high frequency noises
sigf = sgolayfilt(sigj,3,5);

%sigf = sigj;

%% plot filtering result
% figure
% hold on
% plot(sigj)
% plot(sigf)
% hold off
% legend('Raw data','FIR filtering result')

%% EMD component extraction

[imf,residual] = eemd(sigf,0.1,200,1000);


%% calculate the dominant frequencies of different imfs
[m n]=size(imf);

for i = 1:m
   
    ff(i) = mean(instfreq(imf(i,:),fs));
    
end

%% 
Hrlow = 0.75;
Hrhigh = 2;
Rrlow = 0.05;
Rrhigh = 0.75;


HRindex = find(ff>Hrlow & ff<Hrhigh);
RRindex = find(ff>Rrlow & ff<Rrhigh);

%% HR PCA

ftlength = 10000;
df = fs/ftlength;


[coeffH,scoreH,latentH] = pca(imf(min(HRindex):max(HRindex),:)');

hrspectrum = abs(fft(scoreH(:,1),ftlength));

hrspectrumf = sgolayfilt(hrspectrum,3,5);


%%  RR PCA

[coeffR,scoreR,latentR] = pca(imf(min(RRindex):max(RRindex),:)');

rrspectrum = abs(fft(scoreR(:,1),ftlength));

rrspectrumf = sgolayfilt(rrspectrum,3,5);

%% calculate RR HR

hhlow = Hrlow/df;
hhhigh = Hrhigh/df;
rrlow = Rrlow/df;
rrhigh = Rrhigh/df;


HR = (find(hrspectrumf(hhlow:hhhigh) == max(hrspectrumf(hhlow:hhhigh))) +hhlow-2)* df * 60
RR = (find(rrspectrumf(rrlow:rrhigh) == max(rrspectrumf(rrlow:rrhigh))) +rrlow-2)* df * 60

%% plot result
ddf = 0:df:(length(hrspectrum)-1)*df;

figure
hold on
plot(ddf,hrspectrum,'linewidth',2);
plot(ddf,hrspectrumf,'linewidth',2);
plot(ddf,rrspectrum,'linewidth',2);
plot(ddf,rrspectrumf,'linewidth',2);
hold off
xlim([0 2])
box on
xlabel('Frequency (Hz)')
legend('HR','Smooth HR','RR','Smooth RR')
title(['HR = ',num2str(HR),' RR= ',num2str(RR)])