%%
clear
close all
clc

addpath(genpath('code_mfbc'))


%% 

load('fault1_data_n.mat')

% normal_data = zeros(250001,36);
%  for i = 1:12
% %     test = ScopeData.signals(i).values;
%     if size(ScopeData.signals(i).values,2) == 3
%         normal_data(:,(i-1)*3+1:(i-1)*3+3)=ScopeData.signals(i).values;
%     else
%         normal_data(:,(i-1)*3+1)=ScopeData.signals(i).values;
%     end
% % colNames = {'Normal_A','Normal_B','Normal_C'};
% % Normal_Data = array2table(test,'VariableNames',colNames);
%  end
%  
%  
% normal_data(:,[17:18,35:36]) = [];
% 
%  normal_data = downsample(normal_data,100);
 
 
  normal_data = fault1_data_n;
  
 %%
 

 
 
% 
%   figure
%   subplot(121)
%   imagesc( (normal_data))
%   title('raw')
%   subplot(122)
%   imagesc( normc(normal_data));
%   title('normalized')
  
  
  %%
  
  normal_data_ia = abs(hilbert(normal_data));
  bb = 100;
  
  
  figure
  plot(normal_data_ia(bb:end-bb,:));
  title('normal')
  

   %%  leverage score
[U,S,V] = svd(normal_data_ia(bb:end-bb,:),'econ');
lev = vecnorm(U');

figure
plot(lev(bb:end));
title('Leverage Score')


%%

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



 %%
 
 normal_data = downsample(normal_data,100);
 
 
  normal_data_n = normc(normal_data);
  figure
  subplot(121)
  imagesc( (normal_data))
  title('raw')
  subplot(122)
  imagesc( normc(normal_data));
  title('normalized')
  
  
  %%
  
  normal_data_ia = abs(hilbert(normal_data));
  bb = 100;
  
  
  figure
  plot(normal_data_ia(bb:end-bb,:));
  title('normal')
  

   %%  leverage score
[U,S,V] = svd(normal_data_ia(bb:end-bb,:),'econ');
lev = vecnorm(U');

figure
plot(lev(bb:end));
title('Leverage Score')
  
  
  %% Fault 1 
  
opt_anysign = opt_Integerfac_findvert('nonnegative', false, 'affine', false);
r1 = 3;
[That1, Ahat1, status1] = Integerfac_findvert_cpp(normal_data_ia(bb:end-bb,:), r1, [0 1], opt_anysign);

figure; imagesc(That1);title('Binary fault1');
figure; imagesc(Ahat1);title('Coef fault1');
test =  That1 * Ahat1;
figure; plot(test);
figure; plot(normal_data_ia(bb:end-bb,:));


%%

% rest several signals
fault_1_n = normalize(normal_data_ia(bb:end-bb,:),'range');
r1a = 8;
[That1a, Ahat1a, status1a] = Integerfac_findvert_cpp(fault_1_n, r1a, [0 1], opt_anysign);

figure; imagesc(That1a);title('Binary fault1');
figure; imagesc(Ahat1a);title('Coef fault1');
testa =  That1a * Ahat1a;
figure; plot(testa);
figure; plot(fault_1_n);

%% second layer
res_falut1 = fault_1_n - That1a * Ahat1a;
 res_falut1_n = normalize(res_falut1,'range');
r_res1 = 8;
[That1_res1, Ahat1_res1, status1_res1] = Integerfac_findvert_cpp(res_falut1, r_res1, [0 1], opt_anysign);

figure; imagesc(That1_res1);title('Binary fault1');
figure; imagesc(Ahat1_res1);title('Coef fault1');
test1_res1 =  That1_res1 * Ahat1_res1;
figure; plot(test1_res1);
figure; plot(res_falut1);
figure; imagesc(test1_res1);
figure; imagesc(res_falut1);

res2_falut1 = res_falut1 - That1_res1 * Ahat1_res1;

figure; plot(res2_falut1);
