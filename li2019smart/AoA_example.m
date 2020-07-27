clear
close all
clc
warning off

traceN = csvread('Juan_unit17n.csv');
traceE = csvread('Juan_unit17e.csv');

sigN = wdenoise(traceN,'NoiseEstimate','LevelDependent');
sigE = wdenoise(traceE,'NoiseEstimate','LevelDependent');


fs = 500; %sampling rate


total = size(traceN);



%%

wins = 24000;
wine = 24900;

% sigN = sigN(wins:wine);
% sigE = sigE(wins:wine);

% figure;
% subplot(211);
% plot(sigN);
% title 'N';
% xlim([wins wine]);
% ylim([-1.1*abs(max(sigN(wins:wine))) 1.1*max(sigN(wins:wine)) ])
% 
% 
% subplot(212);
% plot(sigE);
% title 'E';
% xlim([wins wine]);
% ylim([-1.1*abs(max(sigE(wins:wine))) 1.1*max(sigE(wins:wine)) ])
 trace = sigN;

%%
 window = 9;
 threshin = std(trace)*(1);
 threshout = std(trace)*(0.5); %1000
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


 
 for i = window + 1:total
   % disp('Entra')
    pos = floor(i -((window/2)));
    
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
        end
    end
 
 end

%%
 
figure

subplot(211);
plot(sigN);
title 'N channel';
xlim([wins wine]);
ylim([-1.1*abs(max(sigN(wins:wine))) 1.1*max(sigN(wins:wine)) ])
hold on
for i = 1:events
    plot([allBeginE(i),allBeginE(i)],[-max(trace),max(trace)],'b');
    plot([allEndE(i),allEndE(i)],[-max(trace),max(trace)],'m');
end
hold off;
xticks([24000 24100 24200 24300 24400 24500 24600 24700 24800 24900])
xticklabels({'48','48.2','48.4','48.6','48.8','49','49.2','49.4','49.6','49.8','50'})
xlabel('Time (s)')
ylabel('Amplitude')


subplot(212);
plot(sigE);
title 'E channel';
xlim([wins wine]);
ylim([-1.1*abs(max(sigE(wins:wine))) 1.1*max(sigE(wins:wine)) ])
hold on
for i = 1:events
    plot([allBeginE(i),allBeginE(i)],[-max(trace),max(trace)],'b');
    plot([allEndE(i),allEndE(i)],[-max(trace),max(trace)],'m');
end
hold off;
xticks([24000 24100 24200 24300 24400 24500 24600 24700 24800 24900])
xticklabels({'48','48.2','48.4','48.6','48.8','49','49.2','49.4','49.6','49.8','50'})
xlabel('Time (s)')
ylabel('Amplitude')

%% SVD part
xx = allBeginE(16:18)+70;
yy = allEndE(16:18)-20;

for ii =1:3
    
    [u,s,v]=svd([sigE(xx(ii):yy(ii)) sigN(xx(ii):yy(ii))]);

scale = 1.5*max(max(trace));
% x=[-coeff(1,1) coeff(1,1)]*scale;
% y=[-coeff(2,1) coeff(2,1)]*scale;
% z=[-coeff(3,1) coeff(3,1)]*scale;

%%

sx=[-v(1,1) v(1,1)]*scale;
sy=[-v(2,1) v(2,1)]*scale;
% sz=[-v(3,1) v(3,1)]*scale;

figure
hold on
scatter(sigN(xx(ii):yy(ii)),sigE(xx(ii):yy(ii)),'k');
xlabel 'N component'
ylabel 'E component'
% plot3(x,y,z,'r','linewidth',2);
plot(sy,sx,'c','linewidth',2);
hold off
title('Angle between N and E directions')
xlim([-1000 1000])
ylim([-1000 1000])
legend('vibration samples','angle vector')
box on
    
end



