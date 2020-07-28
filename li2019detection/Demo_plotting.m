%%

load That0a
load That1a
load That2a
load That3a
load That4a
load That5a
load That6a
load That7a


load That0_res
load That1_res
load That2_res
load That3_res
load That4_res
load That5_res
load That6_res
load That7_res


%%
warning off
aa = 900;
bb = 1080;

That0a = That0a(aa:bb,:);
That1a = That1a(aa:bb,:);
That2a = That2a(aa:bb,:);
That3a = That3a(aa:bb,:);
That4a = That4a(aa:bb,:);
That5a = That5a(aa:bb,:);
That6a = That6a(aa:bb,:);
That7a = That7a(aa:bb,:);


% xx0 = tsne(That0a,'Distance','cityblock');
% xx1 = tsne(That1a,'Distance','cityblock');
% xx2 = tsne(That2a,'Distance','cityblock');
% xx3 = tsne(That3a,'Distance','cityblock');
% xx4 = tsne(That4a,'Distance','cityblock');
% xx5 = tsne(That5a,'Distance','cityblock');
% xx6 = tsne(That6a,'Distance','cityblock');
% xx7 = tsne(That7a,'Distance','cityblock');


% xx0 = tsne(That0a,'Distance','euclidean');
% xx1 = tsne(That1a,'Distance','euclidean');
% xx2 = tsne(That2a,'Distance','euclidean');
% xx3 = tsne(That3a,'Distance','euclidean');
% xx4 = tsne(That4a,'Distance','euclidean');
% xx5 = tsne(That5a,'Distance','euclidean');
% xx6 = tsne(That6a,'Distance','euclidean');
% xx7 = tsne(That7a,'Distance','euclidean');

xx0 = tsne(That0a,'Distance','spearman');
xx1 = tsne(That1a,'Distance','spearman');
xx2 = tsne(That2a,'Distance','spearman');
xx3 = tsne(That3a,'Distance','spearman');
xx4 = tsne(That4a,'Distance','spearman');
xx5 = tsne(That5a,'Distance','spearman');
xx6 = tsne(That6a,'Distance','spearman');
xx7 = tsne(That7a,'Distance','spearman');


%%
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


%%

warning off
aa = 900;
bb = 1080;



%That0_res = That0_res(aa:bb,:);
That1_res = That1_res(aa:bb,:);
That2_res = That2_res(aa:bb,:);
That3_res = That3_res(aa:bb,:);
That4_res = That4_res(aa:bb,:);
That5_res = That5_res(aa:bb,:);
That6_res = That6_res(aa:bb,:);
That7_res = That7_res(aa:bb,:);



% xx0 = tsne(That0a,'Distance','cityblock');
% xx1 = tsne(That1a,'Distance','cityblock');
% xx2 = tsne(That2a,'Distance','cityblock');
% xx3 = tsne(That3a,'Distance','cityblock');
% xx4 = tsne(That4a,'Distance','cityblock');
% xx5 = tsne(That5a,'Distance','cityblock');
% xx6 = tsne(That6a,'Distance','cityblock');
% xx7 = tsne(That7a,'Distance','cityblock');


% xx0 = tsne(That0a,'Distance','euclidean');
% xx1 = tsne(That1a,'Distance','euclidean');
% xx2 = tsne(That2a,'Distance','euclidean');
% xx3 = tsne(That3a,'Distance','euclidean');
% xx4 = tsne(That4a,'Distance','euclidean');
% xx5 = tsne(That5a,'Distance','euclidean');
% xx6 = tsne(That6a,'Distance','euclidean');
% xx7 = tsne(That7a,'Distance','euclidean');

%%

%yy0 = tsne(That0_res);
yy1 = tsne(That1_res);
yy2 = tsne(That2_res);
yy3 = tsne(That3_res);
yy4 = tsne(That4_res);
yy5 = tsne(That5_res);
yy6 = tsne(That6_res);
yy7 = tsne(That7_res);


%%
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

