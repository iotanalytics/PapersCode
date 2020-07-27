close all ; clear all ;

[header, recorddata] = edfread('r01.edf') ;
x0 = recorddata(5,:) ;


x = zeros(size(x0)) ;
for ii = 1:length(x)
    idx = max(1, ii-35):min(length(x), ii+35) ;
    x(ii) = x0(ii) - median(x0(idx)) ;
end
x0 = x ;


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

tfrsq = zeros(500, 1483) ;
tfr = zeros(500, 1483) ;
rtfr = zeros(500, 1483) ;
tfrr = zeros(500, 1483) ;

for qqq = 1: DS

	x = x0(qqq:DS:end) ; %resample(x,1,DS);

	[tfrtmp, ceps, tceps, tfrrtmp, rtfrtmp, tfrsqtmp, tfrtic] = CFPH(x, basicTF, advTF, cepR, P);
	tfr = tfr + tfrtmp ;
	tfrsq = tfrsq + tfrsqtmp ;
	rtfr = rtfr + rtfrtmp ;
	tfrr = tfrr + tfrrtmp ;
	
end

x = resample(x0,1,DS);
[tfrQ, ceps, tceps, tfrrQ, rtfrQ, tfrsqQ, tfrtic] = CFPH(x, basicTF, advTF, cepR, P);

figure(1)
subplot(222)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfrQ, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); axis([-inf inf 0 4]) ;

subplot(221)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.995); axis xy;
xlabel('time (s)'); ylabel('freqnency (Hz)'); axis([-inf inf 0 4]) ;

subplot(223)
imageSQ(0:time_stamp:time_stamp*(size(rtfrtmp,2)-1), tfrtic*basicTF.fs, rtfrtmp, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');  axis([-inf inf 0 4]) ;

subplot(224)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.999); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');  axis([-inf inf 0 4]) ;



