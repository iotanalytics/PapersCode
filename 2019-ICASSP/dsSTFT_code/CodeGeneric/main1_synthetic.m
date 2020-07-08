clear all ; close all ;
initstate(1) ;
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


ff = abs(cumsum(randn(size(x1)))) ;
IF3 = ff./max(abs(ff)) + 4 ;
IF3 = smooth(IF3, 10000) ;
IF3(30001:end) = IF3(30001:end) - t(1:70001)/40 ;
phi = cumsum(IF3) ./ Hz ;
AM3 = smooth(abs(cumsum(randn(size(x1)))./Hz) + 1, 20000) ;
AM3 = AM3 ./ max(AM3) + .9 ;
gg = mod(phi,1);
[a,b] = findpeaks(gg);
s3 = zeros(size(phi)) ;
for ii = 1: length(b)
    idx = max(1,b(ii)-30) : min(b(ii)+30, length(s3)) ;
    s3(idx) = 1 ;
end
x3 = AM3 .* s3 ;
x3(1:30000) = 0 ;
IF3(1:30000) = nan ;


sigma = 0.5 ;%sqrt( var(clean)*10.^( -snrdb /10 ) );
noise = random('T',4,length(x3),1) ;
noise = sigma * noise ; 
var(noise)
snrdb = 20 * log10(std(x1+x2+x3)./std(noise)) ;
fprintf(['snrdb = ',num2str(snrdb),'\n']) ;


x = x1 + x2 + x3 ;
Y = x + noise ;

%noise = randn(size(x)) * 0.5 ;
%snrdb = 20*log(std(x)/std(noise))
x = x(1:20:end) ;
Y = Y(1:20:end) ;
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

[tfrC, cepsC, tcepsC, tfrrC, rtfrC, tfrsqC, tfrtic] = CFPH(x-mean(x), basicTF, advTF, cepR, P);
[tfr, ceps, tceps, tfrr, rtfr, tfrsq, tfrtic] = CFPH(Y-mean(Y), basicTF, advTF, cepR, P);

time_stamp = basicTF.hop/basicTF.fs;




figure ;
subplot(411) ;
plot(t,x3(1:20:end),'k') ; hold on ; set(gca,'fontsize', 18) ; axis([20 80 -inf inf]) ;
set(gca,'fontsize', 18) ; axis([25 65 1 inf]) ;
subplot(412) ;
plot(t, x, 'k') ; axis([25 65 -inf inf]) ;
set(gca,'fontsize', 18) ; 
subplot(413) ;
plot(t,noise(1:20:end),'k') ; axis([25 65 -inf inf]) ;
set(gca,'fontsize', 18) ; 
subplot(414) ;
plot(t,Y,'k') ; axis([25 65 -inf inf]) ;
set(gca,'fontsize', 18) ; xlabel('Time (s)') ;






figure ;
subplot(231)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrC, 0.995); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)'); 

subplot(232)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrsqC, 0.995); axis xy;
xlabel('time (s)'); ylabel('freqnency (s)');

subplot(233)

imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), (1:size(ceps,1))./basicTF.fs, abs(cepsC), .999); axis xy; axis([-inf inf 0 5]);
xlabel('time (s)'); ylabel('quefrency (sec)'); 

subplot(234)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tcepsC, 0.995); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');

subplot(235)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfrrC, 0.998); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');

subplot(236)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfrC, 0.995); axis xy; colormap(1-gray);
%hold on ; plot(t, [ones(3001,1); nan(2000,1)], 'r') ;
xlabel('time (s)'); ylabel('frequency (Hz)');


scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/2]) ;

subplot(131)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, tfr, 0.995); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');

subplot(132)
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfr, 0.995); axis xy; colormap(1-gray);
xlabel('time (s)'); ylabel('frequency (Hz)');


subplot(133) ; hold off;
imageSQ(0:time_stamp:time_stamp*(size(rtfr,2)-1), tfrtic*basicTF.fs, rtfrC, 0.995); axis xy; colormap(1-gray);
hold on ; plot(t, [ones(3001,1); nan(2000,1)], 'r') ;
plot(t, IF2(1:20:end), 'r') ;
plot(t, IF3(1:20:end), 'r') ;
xlabel('time (s)'); ylabel('frequency (Hz)');
