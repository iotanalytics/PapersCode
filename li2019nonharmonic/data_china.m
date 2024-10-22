clear
close all
clc

%%

% load ppg1
% load ppg2
% 
% figure
% subplot(211)
% plot(ppg1)
% subplot(212)
% plot(ppg2)


%%

%sig1 = csvread('TestPPGdata/TestPPGdata1.csv');

sig2 = csvread('TestPPGdata/TestPPGdata2.csv');

%sig3 = csvread('TestPPGdata/TestPPGdata3.csv');

%sig4 = csvread('TestPPGdata/TestPPGdata4.csv');

%%

% figure
% subplot(211)
% plot(sig1(:,1))
% title 'data1 raw'
% xlim([0 length(sig1)])
% subplot(212)
% plot(sig1(:,2))
% title 'data1 filtering'
% xlim([0 length(sig1)])
% 
% 
% figure
% subplot(211)
% plot(sig2(:,1))
% title 'data2 raw'
% xlim([0 length(sig2)])
% subplot(212)
% plot(sig2(:,3))
% title 'data2 filtering'
% xlim([0 length(sig2)])
% 
% 
% figure
% subplot(211)
% plot(sig3(:,1))
% title 'data3 raw'
% xlim([0 length(sig3)])
% subplot(212)
% plot(sig3(:,3))
% title 'data3 filtering'
% xlim([0 length(sig3)])
% 
% 
% 
% figure
% subplot(211)
% plot(sig4(:,1))
% title 'data4 raw'
% xlim([0 length(sig4)])
% subplot(212)
% plot(sig4(:,3))
% title 'data4 filtering'
% xlim([0 length(sig4)])


%%

testppg = sig2(985550:989000,1);
% figure
% plot(testppg)


    N = 2^11;
    x = (0:N-1)/N;
    sig = testppg(1:N)';
    
    numGroup = 3;
    opt.eps = 1e-3;
    opt.res = 0.1;
    opt.freq_range = [0 N/16];
    opt.NG = N;
    opt.dt = 1/N;
    opt.t_sc = 0.5;
    opt.NM = 0;
%     opt.st = round([1 10 ]/16/opt.res)*N/2^11;
%     opt.ed = round([10 20 ]/16/opt.res)*N/2^11;
    opt.st = [5 20 90];
    opt.ed = [100 120 150];


    opt.num_select = numGroup;
    opt.red = 5;
    opt.C = 100;
    opt.rad = 1.5;
    opt.show = 0;
    [insFreq,insAmp,insPhase,comp_select] = insInfo(sig,opt);
    % find peaks
    peaks = zeros(size(comp_select));
    for cnt = 1:numGroup
        peaks(cnt,:) = peakDetection(comp_select(cnt,:),insFreq(cnt,:));
    end
    
    % correct phases
    insPhase = phaseShift(insPhase,peaks);
    
    %%
    opt.maxiter = 100;
    opt.eps_error = 1e-6;
    opt.show = 0;
    opt.iterStyle = 'GS';
    opt.shapeMethod = 1;
    opt.eps_diff = 1e-6;
    opt.ampErrBandWidth = 0;
    opt.numSweep = 1;
    
    switch opt.shapeMethod
        case 1
            opt.para.Ls=1000;
            opt.para.bandWidth = 10;
            opt.para.diffeoMethod = 'nufft';
        case 2
            opt.para.nknots = 20;
            opt.para.knotremoval_factor= 1.0001;
            opt.para.order = 3;
            opt.para.Ls = 1000;
    end
    
    % test example: two components
    [shape,comp] = DeCom_MMD(sig,x,numGroup,insAmp,insFreq,insPhase,opt);
    
    %% Modified to regenerate the signal to test if it is actually accurate enough
    Sig = comp{1} + comp{2}+ comp{3};
    
    N1 = 1; N2 = N;
    Ntotal = length(Sig);
    %%
    pic = figure;
    plot((N1:N2)/Ntotal*480,sig(N1:N2),'b'); axis tight; xlabel('time');%ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/RDSA_fig1_comp1_p2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig1_comp1_p2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    %%
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp{1}(N1:N2),'b'); axis tight; xlabel('time');%ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/RDSA_fig1_comp1_p2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig1_comp1_p2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp{2}(N1:N2),'b'); axis tight; xlabel('time');%ylabel('signal intensity');
   pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/RDSA_fig1_comp2_p2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig1_comp2_p2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp{3}(N1:N2),'b'); axis tight; xlabel('time');%ylabel('signal intensity');
   pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/RDSA_fig1_comp2_p2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig1_comp2_p2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
%     pic = figure;
%     plot((N1:N2)/Ntotal*480,sig(N1:N2)-Sig(N1:N2),'b'); axis tight; xlabel('time');%ylabel('signal intensity');
%     pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
%     saveas(pic,'./results/RDSA_fig1_res_p2.fig');
%     set(gca, 'FontSize', 16);
%     b=get(gca);
%     set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
%     str = './results/RDSA_fig1_res_p2';
%     print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
