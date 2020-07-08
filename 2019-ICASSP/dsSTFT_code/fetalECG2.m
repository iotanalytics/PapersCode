close all ; clear all ;

[header, recorddata] = edfread('r01.edf') ;


DS = 10 ;

basicTF.win = 700; %4096;
basicTF.hop = 20; %441;
basicTF.fs = 1000/DS;
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


time_stamp = basicTF.hop/basicTF.fs;

x0 = recorddata(2,:) ; x = zeros(size(x0)) ;
for ii = 1:length(x)
    idx = max(1, ii-35):min(length(x), ii+35) ;
    x(ii) = x0(ii) - median(x0(idx)) ;
end
x0 = x ; x = resample(x0,1,DS);
[tfrQ2, ceps2, tceps2, tfrrQ2, rtfrQ2, tfrsqQ2, tfrtic2] = CFPH(x, basicTF, advTF, cepR, P);

clear x ; x0 = recorddata(3,:) ; x = zeros(size(x0)) ;
for ii = 1:length(x)
    idx = max(1, ii-35):min(length(x), ii+35) ;
    x(ii) = x0(ii) - median(x0(idx)) ;
end
x0 = x ; x = resample(x0,1,DS);
[tfrQ3, ceps3, tceps3, tfrrQ3, rtfrQ3, tfrsqQ3, tfrtic3] = CFPH(x, basicTF, advTF, cepR, P);


clear x ; x0 = recorddata(4,:) ; x = zeros(size(x0)) ;
for ii = 1:length(x)
    idx = max(1, ii-35):min(length(x), ii+35) ;
    x(ii) = x0(ii) - median(x0(idx)) ;
end
x0 = x ; x = resample(x0,1,DS);
[tfrQ4, ceps4, tceps4, tfrrQ4, rtfrQ4, tfrsqQ4, tfrtic4] = CFPH(x, basicTF, advTF, cepR, P);


clear x ; x0 = recorddata(5,:) ; x = zeros(size(x0)) ;
for ii = 1:length(x)
    idx = max(1, ii-35):min(length(x), ii+35) ;
    x(ii) = x0(ii) - median(x0(idx)) ;
end
x0 = x ; x = resample(x0,1,DS);
[tfrQ5, ceps5, tceps5, tfrrQ5, rtfrQ5, tfrsqQ5, tfrtic5] = CFPH(x, basicTF, advTF, cepR, P);

rmask = (tfrrQ2 + tfrrQ3 + tfrrQ4 + tfrrQ5)./4 ;
maskr = (rtfrQ2 + rtfrQ3 + rtfrQ4 + rtfrQ5)./4 ;


TIME = 0:time_stamp:time_stamp*(size(rtfrQ5,2)-1) ;
FREQ = tfrtic5*basicTF.fs ;

figure(1)
subplot(231); imageSQ(TIME, FREQ, rtfrQ2, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); axis([-inf inf 0 4]) ;
title('abd1 deSST')

subplot(232); imageSQ(TIME, FREQ, rtfrQ3, 0.995); axis xy;
xlabel('time (s)'); ylabel('freqnency (Hz)'); axis([-inf inf 0 4]) ;
title('abd2 deSST')

subplot(233); imageSQ(TIME, FREQ, rtfrQ4, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');  axis([-inf inf 0 4]) ;
title('abd3 deSST')

subplot(234); imageSQ(TIME, FREQ, rtfrQ5, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');  axis([-inf inf 0 4]) ;
title('abd4 deSST')

subplot(235); imageSQ(TIME, FREQ, rmask, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');  axis([-inf inf 0 4]) ;
title('averaged dsSTFT')

subplot(236); imageSQ(TIME, FREQ, maskr, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');  axis([-inf inf 0 4]) ;
title('averaged dsSST')

[c] = CurveExt_M(rtfrQ2', 1);

for ii = 1:size(rtfrQ2,2)
	rtfrQ2(c(ii)-5:c(ii)+5, ii) = 0 ;
	rtfrQ3(c(ii)-5:c(ii)+5, ii) = 0 ;
	rtfrQ4(c(ii)-5:c(ii)+5, ii) = 0 ;
	rtfrQ5(c(ii)-5:c(ii)+5, ii) = 0 ;
end

rtfrQQ = (rtfrQ2 + rtfrQ3 + rtfrQ4 + rtfrQ5) ./ 4 ;

figure(1) ;
subplot(231); imageSQ(TIME, FREQ, rtfrQ2, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); axis([-inf inf 0 4]) ;
title('abd1 deSST')

subplot(232); imageSQ(TIME, FREQ, rtfrQ3, 0.995); axis xy;
xlabel('time (s)'); ylabel('freqnency (Hz)'); axis([-inf inf 0 4]) ;
title('abd2 deSST')

subplot(233); imageSQ(TIME, FREQ, rtfrQ4, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');  axis([-inf inf 0 4]) ;
title('abd3 deSST')

subplot(234); imageSQ(TIME, FREQ, rtfrQ5, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');  axis([-inf inf 0 4]) ;
title('abd4 deSST')

subplot(235); imageSQ(TIME, FREQ, rtfrQQ, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');  axis([-inf inf 0 4]) ;
title('averaged dsSST')

subplot(236); imageSQ(TIME, FREQ, rtfrQ4./(1+tfrrQ3), 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');  axis([-inf inf 0 4]) ;
title('averaged dsSST')

