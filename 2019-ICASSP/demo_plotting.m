
clear
close all
clc


ppg_syn = load('RRSYNTHdata.mat');


%%
figure
subplot(511)
plot(0:1/500:(length(ppg_syn.data(1).ppg.v)-1)/500,ppg_syn.data(1).ppg.v)
xlabel 'Time (s)'
ylabel 'Amplitude'
ylim([-0.1 1.1])
title 'No Modulation'
subplot(512)
plot(0:1/500:(length(ppg_syn.data(2).ppg.v)-1)/500,ppg_syn.data(2).ppg.v)
xlabel 'Time (s)'
ylabel 'Amplitude'
ylim([-0.2 1.3])
title 'Baseline Wander Modulation (BWM)'
subplot(513)
plot(0:1/500:(length(ppg_syn.data(3).ppg.v)-1)/500,ppg_syn.data(3).ppg.v)
xlabel 'Time (s)'
ylabel 'Amplitude'
ylim([-2.2 4])
title 'Amplitude Modulation (AM)'
subplot(514)
plot(0:1/500:(length(ppg_syn.data(4).ppg.v)-1)/500,ppg_syn.data(4).ppg.v)
xlabel 'Time (s)'
ylabel 'Amplitude'
ylim([0 1.1])
title 'Frequency Modulation (FM)'
subplot(515)
plot(0:1/500:(length(ppg_syn.data(5).ppg.v)-1)/500,ppg_syn.data(5).ppg.v)
xlabel 'Time (s)'
ylabel 'Amplitude'
ylim([-0.2 7])
title 'Baseline Wander, Amplitude, Frequency Modulation (BWAFM)'