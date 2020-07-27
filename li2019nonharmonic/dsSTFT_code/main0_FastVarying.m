close all; clear all;


%initstate(1) ;

Hz = 1000 ;
t = [1:100001]' / Hz ;

ff = abs(cumsum(randn(1e5+1,1))) ; IF2 = 3*ff./(max(abs(ff))/2) + pi/2 ;
IF2 = smooth(IF2, 4000) ;
phi = cumsum(IF2) ./ Hz ;
AM2 = smooth(abs(cumsum(randn(size(ff)))./Hz) + 1, 20000) ;
AM2 = AM2 ./ max(AM2) + .9 ;
gg = mod(phi,1);
[a,b] = findpeaks(gg);
b = [1; b; 2*b(end)-b(end-1)] ;
s2 = zeros(size(phi)) ;
for ii = 1: length(b)-1
    idx = b(ii):b(ii+1) ;
    s2(idx) = (idx-b(ii)) ./ (b(ii+1)-b(ii)+1) ;
end
x2 = s2(1:length(AM2)) ;



x = x2 ;

noise = randn(size(x)) * 0.0 ;
snrdb = 20*log(std(x)/std(noise))
x = x(1:20:end) + noise(1:20:end) ;
t = t(1:20:end) ;
Hz = Hz / 20 ;

basicTF.win = 300; %4096;
basicTF.hop = 10; %441;
basicTF.fs = Hz;
basicTF.fr = 0.02;
basicTF.feat = 'SST11';
advTF.num_tap = 1;%num_tap(cc);
advTF.win_type = 'Gauss'; %{'Gaussian','Thomson','multipeak','SWCE'}; %}
advTF.Smo = 1;
advTF.Rej = 0;
advTF.ths = 1E-9;
advTF.HighFreq = 8/50;
advTF.LowFreq = 0.1/50;
advTF.lpc = 0;
cepR.g = 0.1; % g(cc);%0.06;
cepR.Tc=0;
P.num_s = 1;
P.num_c = 1;

[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(x, basicTF, advTF, cepR, P);
% rtfr = medfilt2(rtfr,[1 3]);
time_stamp = basicTF.hop/basicTF.fs;







figure ;
subplot(311) ;
plot(t,x2(1:20:end),'k') ; hold on ; set(gca,'fontsize', 18) ; axis([25 65 -inf inf]) ;
subplot(312) ;
plot(t,IF2(1:20:end),'k') ; hold on ; 
text(56, IF2((end-1)/2)+0.1, '\phi_2''(t)','fontsize',18) ;
plot(t,AM2(1:20:end),'k--') ; hold on ; 
text(t((end-1)/2), AM2((end-1)/2)-0.4, 'A_2(t)','fontsize',18) ;
set(gca,'fontsize', 18) ; axis([20 80 1 inf]) ;
subplot(313) ;
plot(t,x,'color',[.7 .7 .7]) ; hold on ;
plot(t,x2(1:20:end), 'k') ; axis([25 65 -inf inf]) ;
set(gca,'fontsize', 18) ; xlabel('Time') ;






figure ;
subplot(221)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.998); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); 

subplot(222)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsq, 0.998); axis xy;
xlabel('time (s)'); ylabel('freqnency (s)');

subplot(223)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.998); axis xy; colormap(1-gray); axis([0 100 0 6]) ;
xlabel('time (s)'); ylabel('frequency (Hz)'); 

subplot(224)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.998); axis xy; colormap(1-gray);
hold on ; plot(t, [ones(3001,1); nan(2000,1)], 'r') ;
plot(t, IF2(1:20:end), 'r') ; axis([0 100 0 6]) ;
xlabel('time (s)'); ylabel('frequency (Hz)');

