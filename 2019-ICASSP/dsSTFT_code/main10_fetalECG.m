close all ; clear all ;

[header, recorddata] = edfread('r01.edf') ;
x0 = recorddata(5,:) ;
[ann,type,subtype,chan,num,comments] = rdann('r01.edf','qrs');


x = zeros(size(x0)) ;
for ii = 1:length(x)
    idx = max(1, ii-35):min(length(x), ii+35) ;
    x(ii) = x0(ii) - median(x0(idx)) ;
end

ann0 = ann ;
for ii = 1: length(ann0) ;
	[a,b] = max(x(ann(ii)-5:ann(ii)+5)) ;
	ann(ii) = ann(ii) - 5 + b -1 ;
end

fIHR = interp1((ann(2:end)+ann(1:end-1))/2, 1000./(ann(2:end)-ann(1:end-1)), [1:size(recorddata,2)],'pchip','extrap') ;


DS = 20 ;
x = resample(x,1,DS);

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
hold on ; plot([1:length(fIHR)]/1000, fIHR, 'r') ; axis([-inf inf 0 4]) ;
xlabel('time (s)'); ylabel('frequency (Hz)'); 





scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/4]) ;
plot([1:length(recorddata(4,1:4e4))]/1000, x0(1:4e4),'k')
hold on; plot(ann(1:86)/1000, x0(ann(1:86)),'ro','markersize',8) ;
axis([10 40 -inf inf]) ;
set(gca,'fontsize', 18) ; xlabel('Time (s)') ; ylabel('arbitrary unit') ;

