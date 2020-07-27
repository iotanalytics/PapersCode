%% for loop for construct feature matrix
clear;
close all;
clc;
addpath(genpath('code_mfbc'))

load('normal_data_n.mat')

% number of attack type
n = 7;

for i = 1:n
    eval(['load fault', num2str(i),'_data_n'])
end
%%
% for every fault situation

% get the envelope and the difference feature

for i = 1:n
    eval(['fault_this_data_ia = abs(hilbert(fault',num2str(i),'_data_n));']);
    diff_fault_this_data = diff(fault_this_data_ia,1,1);
   
    ub_temp = ([fault_this_data_ia(:,1:15),fault_this_data_ia(:,17:31)]);
    for j = 1:3:28
        ub(:,j) = ((ub_temp(:,j)-ub_temp(:,j+1)).^2+(ub_temp(:,j)-ub_temp(:,j+2)).^2+(ub_temp(:,j+1)-ub_temp(:,j+2)).^2)/3;
        ub(:,j+1) = (max([ub_temp(:,j),ub_temp(:,j+1),ub_temp(:,j+2)],[],2)-min([ub_temp(:,j),ub_temp(:,j+1),ub_temp(:,j+2)],[],2))./max([ub_temp(:,j),ub_temp(:,j+1),ub_temp(:,j+2)],[],2);
        ub(:,j+2) = (abs(ub_temp(:,j)-ub_temp(:,j+1))+abs(ub_temp(:,j)-ub_temp(:,j+2))+abs(ub_temp(:,j+1)-ub_temp(:,j+2)))/3;     
    end

    featureM = [normc(fault_this_data_ia(2:end,:)),normc(diff_fault_this_data), normc(ub(2:end,:))];
    eval(['figure(',num2str(i),')']);
    imagesc(featureM,[-0.04 0.04])
    set(gca,'xtick',[],'ytick',[])
    box on
    bb = 100;
    ylim([bb 2400])
    
    eval(['featureM_',num2str(i),' = featureM;']);
    eval(['fault',num2str(i),'_data_ia = fault_this_data_ia;']);
end


%% leverage score
close all;
for i = 1:n
    eval(['fault_this_data_ia = fault',num2str(i),'_data_ia;']);
    eval(['featureM = featureM_',num2str(i),';']);
    
    %   raw and feature data with noise
    bb = 100;
    noise = 0.003;
    %   norm situation with noise
    [U1,S1,V1] = svd(fault_this_data_ia(bb:end-bb,:)+noise*randn(size(fault_this_data_ia(bb:end-bb,:))),'econ');
    lev1 = vecnorm(U1');
    %   high-dimensional situation with noise
    [U2,S2,V2] = svd(featureM(bb:end-bb,:)+noise*randn(size(featureM(bb:end-bb,:))),'econ');
    lev2 = vecnorm(U2');

    figure(i)
    subplot(1,2,1)
    hold on
    plot(lev2,'linewidth',1.5);
    plot(lev1,'linewidth',1.5);
    title('raw and feature data with noise')
    xlim([0 2302])
    box on
    set(gca,'xtick',[],'ytick',[])
    legend('Feature based','Raw Data based')
    
    % raw and feature data without noise
    [U3,S3,V3] = svd(fault1_data_ia(bb:end-bb,:),'econ');
    lev3 = vecnorm(U3');

    [U4,S4,V4] = svd(featureM(bb:end-bb,:),'econ');
    lev4 = vecnorm(U4');

    subplot(1,2,2)
    hold on
    plot(lev4(bb:end-bb),'linewidth',1.5);
    plot(lev3(bb:end-bb),'linewidth',1.5);
    title('raw and feature data without noise')
    xlim([0 2102])
    box on
    set(gca,'xtick',[],'ytick',[])
    legend('Feature based','Raw Data based')
end

%% BMF for normal data and 7 fault data
% close all the figure windows
close all;
r1a = 11;
for i = 1:n+1
    if i==1
        fault_this_data_ia = abs(hilbert(normal_data_n));
    else
        eval(['fault_this_data_ia = fault',num2str(i-1),'_data_ia;']);
    end
        
    opt_anysign = opt_Integerfac_findvert('nonnegative', false, 'affine', false);

    fault_this_n = normalize(fault_this_data_ia(bb:end-bb,:),'range');

    [That0a, Ahat0a, status0a] = Integerfac_findvert_cpp(fault_this_n, r1a, [0 1], opt_anysign);
%    first layer BMF
    figure(i);
    subplot(1,2,1)
    imagesc(That0a);
    title('normal');
    ylim([800 1200])

    res_falut0 = fault_this_n - That0a * Ahat0a;
    res_falut0_n = normalize(res_falut0,'range');
    [That0_res, Ahat0_res, status0_res] = Integerfac_findvert_cpp(res_falut0, r1a, [0 1], opt_anysign);

%    second layer BMF
    subplot(1,2,2);
    imagesc(That0_res);
    title('normal');
    ylim([800 1200])

    eval(['That',num2str(i-1),'a = That0a;']);
    eval(['That',num2str(i-1),'_res = That0_res;']);
        
end

    
%% attack diagnosis first layer
close all;    
warning off
aa = 900;
bb = 1080;

for i=1:n+1
    eval(['That_this_a = That',num2str(i-1),'a(aa:bb,:);']);
    xx = tsne(That_this_a, 'Distance','spearman');
    eval(['xx',num2str(i-1),'= xx;']);
end

figure
hold on

scatter(xx0(:,1),xx0(:,2),'filled')
scatter(xx1(:,1),xx1(:,2),'filled')
scatter(xx2(:,1),xx2(:,2),'filled')
scatter(xx3(:,1),xx3(:,2),'filled')
scatter(xx4(:,1),xx4(:,2),'filled')
scatter(xx5(:,1),xx5(:,2),'filled')
scatter(xx6(:,1),xx6(:,2),'filled')
scatter(xx7(:,1),xx7(:,2),'filled','k')
legend('Normal','A1','A2','A3','A4','A5','A6','A7')
box on

%% diagnosis second layer
close all;
warning off
aa = 900;
bb = 1080;

for i=1:n+1
    eval(['That_this_res = That',num2str(i-1),'_res(aa:bb,:);']);
    yy = tsne(That_this_res);
    eval(['yy',num2str(i-1),'= yy;']);
end

figure
hold on

scatter(yy0(:,1),yy0(:,2),'w')
scatter(yy1(:,1),yy1(:,2),'filled')
scatter(yy2(:,1),yy2(:,2),'filled')
scatter(yy3(:,1),yy3(:,2),'filled')
scatter(yy4(:,1),yy4(:,2),'filled')
scatter(yy5(:,1),yy5(:,2),'filled')
scatter(yy6(:,1),yy6(:,2),'filled')
scatter(yy7(:,1),yy7(:,2),'filled','k')
legend('Normal','A1','A2','A3','A4','A5','A6','A7')
box on
           