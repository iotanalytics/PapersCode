load sig.mat
%load sig2.mat
load comp.mat


%%
figure
subplot(411)
plot(sig,'k','linewidth',2)
xlim([ length(sig)/2  length(sig)])
ylim([-12 12])
set(gca,'xtick',[],'ytick',[])

subplot(412)
plot(comp{1,1},'k','linewidth',2)
xlim([ length(sig)/2  length(sig)])
ylim([-12 12])
set(gca,'xtick',[],'ytick',[])

subplot(413)
plot(comp{1,2},'k','linewidth',2)
xlim([length(sig)/2  length(sig)])
ylim([-12 12])
set(gca,'xtick',[],'ytick',[])

subplot(414)
plot(sig-Sig,'k','linewidth',2)
xlim([ length(sig)/2  length(sig)])
ylim([-12 12])
set(gca,'xtick',[],'ytick',[])

%%

fs = 27;

t = 0:1/fs : (length(sig)/fs-1/fs);

figure
subplot(311)
plot(t,abs(fft(sig)),'k','linewidth',1.5)
%set(gca,'xtick',[],'ytick',[])
xlim([0 6])
subplot(312)
plot(t,abs(fft(comp{1,1})),'k','linewidth',1.5)
set(gca,'xtick',[],'ytick',[])
xlim([0 6])
subplot(313)
plot(t,abs(fft(comp{1,2})),'k','linewidth',1.5)
set(gca,'xtick',[],'ytick',[])
xlim([0 6])
