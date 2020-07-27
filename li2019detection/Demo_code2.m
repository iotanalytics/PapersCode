load('fault_1.mat')

normal_data = zeros(250001,36);
 for i = 1:12
%     test = ScopeData.signals(i).values;
    if size(ScopeData.signals(i).values,2) == 3
        normal_data(:,(i-1)*3+1:(i-1)*3+3)=ScopeData.signals(i).values;
    else
        normal_data(:,(i-1)*3+1)=ScopeData.signals(i).values;
    end
% colNames = {'Normal_A','Normal_B','Normal_C'};
% Normal_Data = array2table(test,'VariableNames',colNames);
 end
 
 
normal_data(:,[17:18,35:36]) = [];

%  normal_data = downsample(normal_data,100);
 
 %%
 
 load('fault1_data_n.mat')
 
 %%
 figure
 subplot(121)
 imagesc(normal_data)
 %set(gca,'xtick',[],'ytick',[])
 box on
 subplot(122)
 imagesc(fault1_data_n)
 %set(gca,'xtick',[],'ytick',[])
 box on
 
 %% get hilbert transform and difference operation
 
 fault1_data_ia = abs(hilbert(fault1_data_n));
 %bb = 100;
  

  diff_fault1_data = diff(fault1_data_ia,1,1);
  
  %%
%   only get the AC node
  ub_temp = ([fault1_data_ia(:,1:15),fault1_data_ia(:,17:31)]);
  
  for i = 1:3:28
     
     ub(:,i) = ((ub_temp(:,i)-ub_temp(:,i+1)).^2+(ub_temp(:,i)-ub_temp(:,i+2)).^2+(ub_temp(:,i+1)-ub_temp(:,i+2)).^2)/3;
     ub(:,i+1) = (max([ub_temp(:,i),ub_temp(:,i+1),ub_temp(:,i+2)],[],2)-min([ub_temp(:,i),ub_temp(:,i+1),ub_temp(:,i+2)],[],2))./max([ub_temp(:,i),ub_temp(:,i+1),ub_temp(:,i+2)],[],2);
     ub(:,i+2) = (abs(ub_temp(:,i)-ub_temp(:,i+1))+abs(ub_temp(:,i)-ub_temp(:,i+2))+abs(ub_temp(:,i+1)-ub_temp(:,i+2)))/3;
      
      
  end
  
  %%
%   get the feature matrix featureM, why 2:end ?? aaronli
  
  featureM = [normc(fault1_data_ia(2:end,:)),normc(diff_fault1_data), normc(ub(2:end,:))];
  figure
  imagesc(featureM,[-0.04 0.04])
  set(gca,'xtick',[],'ytick',[])
  box on
  bb = 100;
  ylim([bb 2400])
%   figure(2);plot(normal_data_ia(bb:end-bb,:));
%   title('normal')
  
  
  %%  leverage score - normal data shall have no peaks, while other fault data will have peaks
  %% raw and feature data with noise
  bb = 100;
%   norm situation with 0.003 noise
[U1,S1,V1] = svd(fault1_data_ia(bb:end-bb,:)+0.003*randn(size(fault1_data_ia(bb:end-bb,:))),'econ');
lev1 = vecnorm(U1');
%   high-dimensional situation with 0.003 noise
[U2,S2,V2] = svd(featureM(bb:end-bb,:)+0.003*randn(size(featureM(bb:end-bb,:))),'econ');
lev2 = vecnorm(U2');

figure
hold on
plot(lev2,'linewidth',1.5);
plot(lev1,'linewidth',1.5);
xlim([0 2302])
box on
set(gca,'xtick',[],'ytick',[])
legend('Feature based','Raw Data based')

%% raw and feature data without noise
[U1,S1,V1] = svd(fault1_data_ia(bb:end-bb,:),'econ');
lev1 = vecnorm(U1');

[U2,S2,V2] = svd(featureM(bb:end-bb,:),'econ');
lev2 = vecnorm(U2');

figure
hold on
plot(lev2(bb:end-bb),'linewidth',1.5);
plot(lev1(bb:end-bb),'linewidth',1.5);
xlim([0 2102])
box on
set(gca,'xtick',[],'ytick',[])
legend('Feature based','Raw Data based')



%%
% 
%  %% Fault 1 
%   
% opt_anysign = opt_Integerfac_findvert('nonnegative', false, 'affine', false);
% r1 = 8;
% [That1, Ahat1, status1] = Integerfac_findvert_cpp(fault1_data_ia(bb:end-bb,:), r1, [0 1], opt_anysign);
% 
% figure; imagesc(That1);title('Binary fault1');
% figure; imagesc(Ahat1);title('Coef fault1');
% test =  That1 * Ahat1;
% figure; plot(test);
% %figure; plot(normal_data_ia(bb:end-bb,:));
% 
% 
% %%
% 
% % rest several signals
% fault_1_n = normalize(fault1_data_ia(bb:end-bb,:),'range');
% r1a = 8;
% [That1a, Ahat1a, status1a] = Integerfac_findvert_cpp(fault_1_n, r1a, [0 1], opt_anysign);
% 
% figure; imagesc(That1a);title('Binary fault1');
% figure; imagesc(Ahat1a);title('Coef fault1');
% testa =  That1a * Ahat1a;
% figure; plot(testa);
% %figure; plot(fault_1_n);
% 
% %% second layer
% res_falut1 = fault_1_n - That1a * Ahat1a;
%  res_falut1_n = normalize(res_falut1,'range');
% r_res1 = 8;
% [That1_res1, Ahat1_res1, status1_res1] = Integerfac_findvert_cpp(res_falut1, r_res1, [0 1], opt_anysign);
% 
% figure; imagesc(That1_res1);title('Binary fault1');
% figure; imagesc(Ahat1_res1);title('Coef fault1');
% test1_res1 =  That1_res1 * Ahat1_res1;
% figure; plot(test1_res1);
% figure; plot(res_falut1);
% figure; imagesc(test1_res1);
% figure; imagesc(res_falut1);
% 
% res2_falut1 = res_falut1 - That1_res1 * Ahat1_res1;
% 
% %figure; plot(res2_falut1);
% 
% 
% %% fault2
% 
%  load('fault2_data_n.mat')
%  fault2_data_ia = abs(hilbert(fault2_data_n));
% 
% opt_anysign = opt_Integerfac_findvert('nonnegative', false, 'affine', false);
% r1 = 8;
% [That1, Ahat1, status1] = Integerfac_findvert_cpp(fault2_data_ia(bb:end-bb,:), r1, [0 1], opt_anysign);
% 
% figure; imagesc(That1);title('Binary fault1');
% figure; imagesc(Ahat1);title('Coef fault1');
% test =  That1 * Ahat1;
% figure; plot(test);
% %figure; plot(normal_data_ia(bb:end-bb,:));
% 
% 
% %%
% 
% % rest several signals
% fault_2_n = normalize(fault2_data_ia(bb:end-bb,:),'range');
% r1a = 8;
% [That2a, Ahat2a, status2a] = Integerfac_findvert_cpp(fault_2_n, r1a, [0 1], opt_anysign);
% 
% % figure; imagesc(That2a);title('Binary fault1');
% figure; imagesc(Ahat2a);title('Coef fault1');
% % test2a =  That2a * Ahat2a;
% % figure; plot(test2a);
% %figure; plot(fault_1_n);
% 
% %% second layer
% res_falut2 = fault_2_n - That2a * Ahat2a;
%  res_falut2_n = normalize(res_falut2,'range');
% r_res2 = 8;
% [That2_res1, Ahat2_res1, status2_res1] = Integerfac_findvert_cpp(res_falut2, r_res2, [0 1], opt_anysign);
% 
% figure; imagesc(That2_res1);title('Binary fault1');
% figure; imagesc(Ahat2_res1);title('Coef fault1');
% test2_res1 =  That2_res1 * Ahat2_res1;
% figure; plot(test2_res1);
% figure; plot(res_falut2);
% figure; imagesc(test2_res1);
% figure; imagesc(res_falut2);
% 
% res2_falut1 = res_falut2 - That2_res1 * Ahat2_res1;
% 
% %figure; plot(res2_falut1);
% 
% 
% %%
% 
% %% fault2
% 
%  load('fault3_data_n.mat')
%  fault3_data_ia = abs(hilbert(fault3_data_n));
% 
% opt_anysign = opt_Integerfac_findvert('nonnegative', false, 'affine', false);
% r1 = 8;
% [That3, Ahat3, status3] = Integerfac_findvert_cpp(fault3_data_ia(bb:end-bb,:), r1, [0 1], opt_anysign);
% 
% figure; imagesc(That3);title('Binary fault1');
% figure; imagesc(Ahat3);title('Coef fault1');
% test =  That3 * Ahat3;
% %figure; plot(test);
% %figure; plot(normal_data_ia(bb:end-bb,:));
% 
% 
% %%
% 
% % rest several signals
% fault_3_n = normalize(fault3_data_ia(bb:end-bb,:),'range');
% r1a = 16;
% [That3a, Ahat3a, status3a] = Integerfac_findvert_cpp(fault_3_n, r1a, [0 1], opt_anysign);
% 
% % figure; imagesc(That2a);title('Binary fault1');
% figure; imagesc(Ahat3a);title('Coef fault1');
% % test2a =  That2a * Ahat2a;
% % figure; plot(test2a);
% %figure; plot(fault_1_n);
% 
% %% second layer
% res_falut3 = fault_3_n - That3a * Ahat3a;
%  res_falut3_n = normalize(res_falut3,'range');
% r_res2 = 16;
% [That3_res1, Ahat3_res1, status3_res1] = Integerfac_findvert_cpp(res_falut3, r_res2, [0 1], opt_anysign);
% 
% %figure; imagesc(That3_res1);title('Binary fault1');
% figure; imagesc(Ahat3_res1);title('Coef fault1');
% %test3_res1 =  That3_res1 * Ahat3_res1;
% % figure; plot(test2_res1);
% % figure; plot(res_falut2);
% % figure; imagesc(test2_res1);
% % figure; imagesc(res_falut2);
% 
% %res3_falut1 = res_falut3 - That3_res1 * Ahat3_res1;
% 
% %figure; plot(res2_falut1);
% 
% 

%%

r1a = 11;

load('normal_data_n.mat')
fault0_data_ia = abs(hilbert(normal_data_n));

opt_anysign = opt_Integerfac_findvert('nonnegative', false, 'affine', false);

fault_0_n = normalize(fault0_data_ia(bb:end-bb,:),'range');

[That0a, Ahat0a, status0a] = Integerfac_findvert_cpp(fault_0_n, r1a, [0 1], opt_anysign);

figure(1);
imagesc(That0a);
title('normal');
ylim([800 1200])

res_falut0 = fault_0_n - That0a * Ahat0a;
res_falut0_n = normalize(res_falut0,'range');
[That0_res, Ahat0_res, status0_res] = Integerfac_findvert_cpp(res_falut0, r1a, [0 1], opt_anysign);



% figure(1);
% imagesc(That0_res);
% title('normal');
% ylim([800 1200])


load('fault1_data_n.mat')
fault1_data_ia = abs(hilbert(fault1_data_n));

opt_anysign = opt_Integerfac_findvert('nonnegative', false, 'affine', false);

fault_1_n = normalize(fault1_data_ia(bb:end-bb,:),'range');

[That1a, Ahat1a, status1a] = Integerfac_findvert_cpp(fault_1_n, r1a, [0 1], opt_anysign);

figure(2);
imagesc(That1a);
title('normal');
ylim([800 1200])

res_falut1 = fault_1_n - That1a * Ahat1a;
res_falut1_n = normalize(res_falut1,'range');
[That1_res, Ahat1_res, status1_res] = Integerfac_findvert_cpp(res_falut1, r1a, [0 1], opt_anysign);

% figure(2)
% imagesc(That1_res);
% title('fault1');
% ylim([800 1200])


load('fault2_data_n.mat')
fault2_data_ia = abs(hilbert(fault2_data_n));

opt_anysign = opt_Integerfac_findvert('nonnegative', false, 'affine', false);

fault_2_n = normalize(fault2_data_ia(bb:end-bb,:),'range');

[That2a, Ahat2a, status2a] = Integerfac_findvert_cpp(fault_2_n, r1a, [0 1], opt_anysign);

figure(3);
imagesc(That0a);
title('fault2');
ylim([800 1200])

res_falut2 = fault_2_n - That2a * Ahat2a;
res_falut2_n = normalize(res_falut2,'range');
[That2_res, Ahat2_res, status2_res] = Integerfac_findvert_cpp(res_falut2, r1a, [0 1], opt_anysign);

% figure(3)
% imagesc(That2_res);
% title('fault2');
% ylim([800 1200])


load('fault3_data_n.mat')
fault3_data_ia = abs(hilbert(fault3_data_n));

opt_anysign = opt_Integerfac_findvert('nonnegative', false, 'affine', false);

fault_3_n = normalize(fault3_data_ia(bb:end-bb,:),'range');

[That3a, Ahat3a, status3a] = Integerfac_findvert_cpp(fault_3_n, r1a, [0 1], opt_anysign);

figure(1);
imagesc(That0a);
title('normal');
ylim([800 1200])

res_falut3 = fault_3_n - That3a * Ahat3a;
res_falut3_n = normalize(res_falut3,'range');
[That3_res, Ahat3_res, status3_res] = Integerfac_findvert_cpp(res_falut3, r1a, [0 1], opt_anysign);

figure(4)
imagesc(That3_res);
title('fault3');
ylim([800 1200])


load('fault4_data_n.mat')
fault4_data_ia = abs(hilbert(fault4_data_n));

opt_anysign = opt_Integerfac_findvert('nonnegative', false, 'affine', false);

fault_4_n = normalize(fault4_data_ia(bb:end-bb,:),'range');

[That4a, Ahat4a, status4a] = Integerfac_findvert_cpp(fault_4_n, r1a, [0 1], opt_anysign);

res_falut4 = fault_4_n - That4a * Ahat4a;
res_falut4_n = normalize(res_falut4,'range');
[That4_res, Ahat4_res, status4_res] = Integerfac_findvert_cpp(res_falut4, r1a, [0 1], opt_anysign);

figure(5)
imagesc(That4_res);
title('fault4');
ylim([800 1200])


load('fault5_data_n.mat')
fault5_data_ia = abs(hilbert(fault5_data_n));

opt_anysign = opt_Integerfac_findvert('nonnegative', false, 'affine', false);

fault_5_n = normalize(fault5_data_ia(bb:end-bb,:),'range');

[That5a, Ahat5a, status5a] = Integerfac_findvert_cpp(fault_5_n, r1a, [0 1], opt_anysign);

res_falut5 = fault_5_n - That5a * Ahat5a;
res_falut5_n = normalize(res_falut5,'range');
[That5_res, Ahat5_res, status5_res] = Integerfac_findvert_cpp(res_falut5, r1a, [0 1], opt_anysign);
% 
figure(6)
imagesc(That5_res);
title('fault5');
ylim([800 1200])


load('fault6_data_n.mat')
fault6_data_ia = abs(hilbert(fault6_data_n));

opt_anysign = opt_Integerfac_findvert('nonnegative', false, 'affine', false);

fault_6_n = normalize(fault6_data_ia(bb:end-bb,:),'range');
[That6a, Ahat6a, status6a] = Integerfac_findvert_cpp(fault_6_n, r1a, [0 1], opt_anysign);

res_falut6 = fault_6_n - That6a * Ahat6a;
res_falut6_n = normalize(res_falut6,'range');
[That6_res, Ahat6_res, status6_res] = Integerfac_findvert_cpp(res_falut6, r1a, [0 1], opt_anysign);

figure(7)
imagesc(That6_res);
title('fault6');
ylim([800 1200])


load('fault7_data_n.mat')
fault7_data_ia = abs(hilbert(fault7_data_n));

opt_anysign = opt_Integerfac_findvert('nonnegative', false, 'affine', false);

fault_7_n = normalize(fault7_data_ia(bb:end-bb,:),'range');

[That7a, Ahat7a, status7a] = Integerfac_findvert_cpp(fault_7_n, r1a, [0 1], opt_anysign);

res_falut7 = fault_7_n - That7a * Ahat7a;
res_falut7_n = normalize(res_falut7,'range');
[That7_res, Ahat7_res, status7_res] = Integerfac_findvert_cpp(res_falut7, r1a, [0 1], opt_anysign);

figure(8)
imagesc(That7_res);
title('fault7');
ylim([800 1200])

%%
figure
subplot(421)
imagesc(That0_res);
ylim([800 1100])
set(gca,'xtick',[],'ytick',[])
colormap 'bone'
subplot(422)
imagesc(That1_res);
ylim([800 1100])
set(gca,'xtick',[],'ytick',[])

subplot(423)
imagesc(That2_res);
ylim([800 1200])
set(gca,'xtick',[],'ytick',[])
colormap 'bone'
subplot(424)
imagesc(That3_res);
ylim([800 1100])
set(gca,'xtick',[],'ytick',[])

subplot(425)
imagesc(That4_res);
ylim([800 1200])
set(gca,'xtick',[],'ytick',[])
colormap 'bone'
subplot(426)
imagesc(That5_res);
ylim([800 1100])
set(gca,'xtick',[],'ytick',[])

subplot(427)
imagesc(That6_res);
ylim([800 1100])
set(gca,'xtick',[],'ytick',[])
colormap 'bone'
subplot(428)
imagesc(That7_res);
ylim([800 1100])
set(gca,'xtick',[],'ytick',[])

%%
figure
subplot(421)
imagesc(That0a);
ylim([800 1100])
set(gca,'xtick',[],'ytick',[])
colormap 'bone'
subplot(422)
imagesc(That1a);
ylim([800 1100])
set(gca,'xtick',[],'ytick',[])

subplot(423)
imagesc(That2a);
ylim([800 1200])
set(gca,'xtick',[],'ytick',[])
colormap 'bone'
subplot(424)
imagesc(That3a);
ylim([800 1100])
set(gca,'xtick',[],'ytick',[])

subplot(425)
imagesc(That4a);
ylim([800 1200])
set(gca,'xtick',[],'ytick',[])
colormap 'bone'
subplot(426)
imagesc(That5a);
ylim([800 1100])
set(gca,'xtick',[],'ytick',[])

subplot(427)
imagesc(That6a);
ylim([800 1100])
set(gca,'xtick',[],'ytick',[])
colormap 'bone'
subplot(428)
imagesc(That7a);
ylim([800 1100])
set(gca,'xtick',[],'ytick',[])









