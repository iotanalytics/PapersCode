clear
close all
clc
warning off
% LAB_61006_E-20180429_190343-20180429_190356.csv 3
% LAB_61006_N-20180429_185259-20180429_185316.csv 2
% LAB_61006_Z-20180429_190829-20180429_190838.csv 1

%traceZ = csvread('LAB_61006_Z-20180429_190829-20180429_190838.csv')
traceZ = csvread('Maria6F2.csv');
traceZ = traceZ';

trace = wdenoise(traceZ,'NoiseEstimate','LevelDependent');

fs = 500; %sampling rate

wins = 0;
wine = 48600;


total = size(traceZ)

%%
 window = 6;
 threshin = std(trace)*(4);
 threshout = std(trace)*(2); %1000
 plot(trace)
 %hold on
 xlim([wins wine]);
 sw = 0;
 countout=0;
 beginE = 0;%zeros(100);
 endE   = 0;%zeros(100);
 
 allBeginE = zeros(100);
 allEndE   = zeros(100);
 events = 0;
 fileID  = fopen('mariaTrain33.txt','w');
 fileID2 = fopen('mariaTest33.txt','w');
 for i = window + 1:total
   % disp('Entra')
    pos = floor(i - ((window/2)));
    
    m2(pos) = moment(trace(i - window : i),2);
    if(m2(pos)>threshin && sw ==0)
        beginE = pos;
        sw = 1;
    end
    if( sw==1 && m2(pos)<threshout )
       countout=countout+1;
    else
       countout=0;
    end
    if( countout>5 )
        endE = pos;
        sw = 0;
        countout=0;
        if(endE-beginE>100)
           events=events+1;
           allBeginE(events)= beginE;
           allEndE(events)  = endE;
           signalE = zeros(endE-beginE);
           signalE = trace(beginE:endE);
           [signalNor] = normalizationE(signalE);
           [stds,entropyS,pksre,maxPeak,locMaxPeak,before,after] = timefeatures(signalNor);
           [spec,centroid,spks,slogs] = frequencyfeatures(signalNor);
           wf = randi([1 100]);
           if(wf < 80)
             fprintf(fileID,'%d ',endE-beginE,stds,entropyS,pksre,maxPeak,locMaxPeak,before,after,spec,centroid,spks);
             fprintf(fileID,'%d',slogs);
             fprintf(fileID,'\n');
           else
             fprintf(fileID,'%d ',endE-beginE,stds,entropyS,pksre,maxPeak,locMaxPeak,before,after,spec,centroid,spks);
             fprintf(fileID,'%d',slogs);
             fprintf(fileID2,'\n');
           end
        end
    end
 
 end

fclose(fileID);
fclose(fileID2);

figure;
subplot(211);
plot(traceZ);
title 'Original';
xlim([wins wine]);


subplot(212);
plot(trace);
title 'Filtered';
xlim([wins wine]);
hold on
for i = 1:events
    plot([allBeginE(i),allBeginE(i)],[-max(trace),max(trace)],'b');
    plot([allEndE(i),allEndE(i)],[-max(trace),max(trace)],'m');
end
hold off;

%plot(m2)

%%
% figure;
% %result1
% begin = [4480,4840,5143,5455,5780,6100];
% %result2
% %begin = [4145,4767,4466,5060,5360,5665];
% %result3
% %begin = [2080,2380,2980,2680,3280,3580];
% plotp = [321,322,323,324,325,326,327,328];
% le = 250;
% 
% for c = 1:6
%     step1 = trace(begin(c):begin(c)+le);
%     x = plotp(c);
%     subplot(x)
%     %step1 = fft(step1);
%     plot(step1)
%     if mod(c,2)==0
%        title 'Left'
%     else
%        title 'Right'
%     end
%     xlim([0 le]) 
%     ylim([-10000 10000]) 
% end
% 
% figure;
% 
% for c = 1:6
%     step1 = real(begin(c):begin(c)+le);
%     x = plotp(c);
%     subplot(x)
%     %step1 = fft(step1);
%     plot(step1)
%     if mod(c,2)==0
%        title 'Left'
%     else
%        title 'Right'
%     end
%     
%     xlim([0 le])  
%     ylim([-10000 10000]) 
% end