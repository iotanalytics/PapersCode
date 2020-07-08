close all ; clear all ;

load Fwave
x = x(1:end/6) ;
x0 = x ;
trend = zeros(size(x)) ;
for ii = 1:length(x0)
    idx = max(1, ii-50):min(length(x), ii+50) ;
    trend(ii) = median(x0(idx)) ;
end

x = x0 - trend ;  

x = resample(x,1,2);

basicTF.win = 200; %4096;
basicTF.hop = 20; %441;
basicTF.fs = 100;
basicTF.fr = 0.02;
basicTF.feat = 'SST11';
advTF.num_tap = 1;%num_tap(cc);
advTF.win_type = 'Gauss'; %{'Gaussian','Thomson','multipeak','SWCE'}; %}
advTF.Smo = 1;
advTF.Rej = 0;
advTF.ths = 1E-9;
advTF.HighFreq = 10/100;
advTF.LowFreq = 0.1/100;
advTF.lpc = 0;
cepR.g = 0.3; % g(cc);%0.06;
cepR.Tc=0;
P.num_s = 1;
P.num_c = 1;


[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(x, basicTF, advTF, cepR, P);
time_stamp = basicTF.hop/basicTF.fs;





figure(1)
subplot(221)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.995); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); 

subplot(222)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
xlabel('time (s)'); ylabel('freqnency (Hz)');

subplot(223)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');  axis([-inf inf 0 4]) ;

subplot(224)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); 




