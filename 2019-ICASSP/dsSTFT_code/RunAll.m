scrsz = get(0,'ScreenSize');

main1_synthetic
set(1,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
export_fig('Main1','-transparent') ;

figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)*3/4]) ;
subplot(3,1,1) ;
plot(t,x) ; axis tight ; set(gca,'fontsize', 18) ; xlabel('Time (s)') ;
subplot(3,1,2) ;
plot(t,x1(1:20:end)) ; hold on ; set(gca,'fontsize', 18) ; xlabel('Time (s)') ; axis tight ;
subplot(3,1,3) ;
plot(t,x2(1:20:end)) ; hold on ; set(gca,'fontsize', 18) ; xlabel('Time (s)') ; axis tight ;
export_fig('Main1sig','-transparent') ;
close all

	%%%
main2_respiration
set(1,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
export_fig('Main2','-transparent') ;

figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/4]) ;
plot([1:length(x1)]'/fs, x1, 'k') ; axis tight ; set(gca,'fontsize', 18) ; xlabel('Time (s)') ;
ylabel('arbitrary unit')
export_fig('Main2sig','-transparent') ;
close all


	%%%
main3_ECG
set(1,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
export_fig('Main3','-transparent') ;

figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/4]) ;
plot([1:length(x)]'/basicTF.fs, x, 'k') ; axis tight ; set(gca,'fontsize', 18) ; xlabel('Time (s)') ; axis([0 30 -inf inf]) ;
ylabel('arbitrary unit')
export_fig('Main3sig','-transparent') ;
close all

	%%%
main4_clarinet
set(1,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
export_fig('Main4','-transparent') ;

figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/4]) ;
plot([1:length(x)]'/basicTF.fs, x, 'k') ; axis tight ; set(gca,'fontsize', 18) ; xlabel('Time (s)') ;
ylabel('arbitrary unit')
export_fig('Main4sig','-transparent') ;
close all

	%%%
main5_piano_octave
set(1,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
export_fig('Main5','-transparent') ;

figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/4]) ;
plot([1:length(x)]'/basicTF.fs, x, 'k') ; axis tight ; set(gca,'fontsize', 18) ; xlabel('Time (s)') ;
ylabel('arbitrary unit')
export_fig('Main5sig','-transparent') ;
close all


	%%%
main6_piano_piece
set(1,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
export_fig('Main6','-transparent') ;
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/4]) ;
plot([1:length(x)]'/basicTF.fs, x, 'k') ; axis tight ; set(gca,'fontsize', 18) ; xlabel('Time (s)') ;
ylabel('arbitrary unit')
export_fig('Main6sig','-transparent') ;
close all

	%%%
main7_ppg2comp_heart_breath
set(1,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
export_fig('Main7','-transparent') ;

figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/4]) ;
plot([1:length(signal.pleth.y)]'/basicTF.fs, signal.pleth.y, 'k') ; axis tight ; set(gca,'fontsize', 18) ; xlabel('Time (s)') ; axis([150 250 -inf inf]) ;
ylabel('arbitrary unit')
export_fig('Main7sig','-transparent') ;
close all

	%%%
main8_ppg3comp_motion
set(1,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
export_fig('Main8','-transparent') ;

figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/4]) ;
plot([1:length(sig(2,:))]'/basicTF.fs, sig(2,:), 'k') ; axis tight ; set(gca,'fontsize', 18) ; xlabel('Time (s)') ;
ylabel('arbitrary unit')
export_fig('Main8sig','-transparent') ;
axis([0 100 -inf inf])
export_fig('Main8sigZoom','-transparent') ;
close all

main9_noncontact_ppg
set(1,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]) ;
export_fig('Main9','-transparent') ;

figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/4]) ;
plot([1:size(x,1)]'/basicTF.fs, x(:,1) - mean(x(:,1)), 'k') ; axis tight ; set(gca,'fontsize', 18) ; xlabel('Time (s)') ;
ylabel('arbitrary unit') ; axis([0 inf -inf 0.015]) ;
export_fig('Main9sig','-transparent') ;
close all
