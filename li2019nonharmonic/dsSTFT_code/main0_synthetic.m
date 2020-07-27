close all; clear all;


%initstate(1) ;
dd = zeros(1e5+1,1) ;
dd(1:1000:end) = 1 ;
aa = exp(-([-2e4:8e4]/1e5).^2) ;
[h,dh] = hermf(201,1,6) ;
h = h ./ max(h) ;

Hz = 1000 ;
t = [1:100001]' / Hz ;
x1 = 3*conv(aa'.*dd,h,'same') ; % 1Hz
x1(60001:end) = 0 ;


ff = abs(cumsum(randn(size(x1)))) ; IF2 = ff./(max(abs(ff))/2) + pi/2 ;
IF2 = smooth(IF2, 10000) ;
phi = cumsum(IF2) ./ Hz ;
AM2 = smooth(abs(cumsum(randn(size(x1)))./Hz) + 1, 20000) ;
AM2 = AM2 ./ max(AM2) + .9 ;
gg = mod(phi,1);
[a,b] = findpeaks(gg);
b = [1; b; 2*b(end)-b(end-1)] ;
s2 = zeros(size(phi)) ;
for ii = 1: length(b)-1
    idx = b(ii):b(ii+1) ;
    s2(idx) = (idx-b(ii)) ./ (b(ii+1)-b(ii)+1) ;
end
x2 = AM2 .* s2(1:length(AM2)) ;



%ff = abs(cumsum(randn(size(x1)))) ; IF2 = ff./max(abs(ff)) + pi/2 ;
%IF2 = smooth(IF2, 10000) ;
%phi = cumsum(IF2) ./ Hz ;
%AM2 = smooth(abs(cumsum(randn(size(x1)))./Hz) + 1, 20000) ;
%AM2 = AM2 ./ max(AM2) + .9 ;
%x2 = AM2 .* mod(phi, 1) ;

x = x1 + x2 ;

%%

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
subplot(411) ;
plot(t,x1(1:20:end),'k') ; hold on ; set(gca,'fontsize', 18) ; axis([25 65 -inf inf]) ;
subplot(412) ;
plot(t,x2(1:20:end),'k') ; hold on ; set(gca,'fontsize', 18) ; axis([25 65 -inf inf]) ;
subplot(413) ;
plot(t,IF2(1:20:end),'k') ; hold on ; 
text(56, IF2((end-1)/2)+0.1, '\phi_2''(t)','fontsize',18) ;
plot(t,AM2(1:20:end),'k--') ; hold on ; 
text(t((end-1)/2), AM2((end-1)/2)-0.4, 'A_2(t)','fontsize',18) ;
set(gca,'fontsize', 18) ; axis([20 80 1 inf]) ;
subplot(414) ;
plot(t,x,'color',[.7 .7 .7]) ; hold on ;
plot(t, x1(1:20:end) + x2(1:20:end), 'k') ; axis([25 65 -inf inf]) ;
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

