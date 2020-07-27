initstate(11)


    %% setup data
dd = zeros(1e5+1,1) ;
dd(1:1000:end) = 1 ;
aa = exp(-([-2e4:8e4]/1e5).^2) ;
[h,dh] = hermf(201,1,6) ;
h = h ./ max(h) ;

Hz = 1000 ;
t = [1:100001]' / Hz ;
x1 = 3*conv(aa'.*dd,h,'same') ; % 1Hz



ff = abs(cumsum(randn(size(x1)))) ; IF2 = ff./(max(abs(ff))/2) + pi/2 ;
IF2 = smooth(IF2, 10000) ;
phi = cumsum(IF2) ./ Hz ;
AM2 = smooth(abs(cumsum(randn(size(x1)))./Hz) + 1, 20000) ;
AM2 = AM2 ./ max(AM2) + .9 ;
gg = mod(phi,1);
[a,b] = findpeaks(gg);
b = [1; b; 2*b(end)-b(end-1)] ;
s2 = zeros(size(phi)) ;
for ii = 1: length(b)-1
    idx = b(ii):b(ii+1) ;
    s2(idx) = (idx-b(ii)) ./ (b(ii+1)-b(ii)+1) ;
end
x2 = AM2 .* s2(1:length(AM2)) ;



sigma = 0 ;
noise = random('T',4,length(x1),1) ;
noise = sigma * noise ;
var(noise)
snrdb = 20 * log10(std(x1+x2)./std(noise)) ;
fprintf(['snrdb = ',num2str(snrdb),'\n']) ;


x = x1 + x2 ;
Y = x + noise ;

%noise = randn(size(x)) * 0.5 ;
%snrdb = 20*log(std(x)/std(noise))
x = x(1:20:end) ;
Y = Y(1:20:end) ;
t = t(1:20:end) ;
Hz = Hz / 20 ;
time_stamp = basicTF.hop/basicTF.fs;




%===================================
	% run de-shape SST
basicTF.win = 301; 
basicTF.hop = 10; 
basicTF.fs = Hz;
basicTF.fr = 0.02;
basicTF.feat = 'SST11';
advTF.num_tap = 1;
advTF.win_type = 'Gauss'; 
advTF.Smo = 1;
advTF.Rej = 0;
advTF.ths = 1E-9;
advTF.HighFreq = 8/50;
advTF.LowFreq = 0.1/50;
advTF.lpc = 0;
cepR.g = 0.1; % g(cc);%0.06;
cepR.Tc=0;
P.num_s = 1;
P.num_c = 1;

[h,dh] = hermf(basicTF.win, 1, 6) ;
[tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFTbase(x-mean(x), 0, 0.5, 0.0005, 1, h(1,:)', dh(1,:)', 0, 0) ;
[tfrC, cepsC, tcepsC, tfrrC, rtfrC, tfrsqC, tfrticC] = CFPH(x-mean(x), basicTF, advTF, cepR, P);


tt = t(1:basicTF.hop:end) ;
tt = tt(1:size(rtfrC,2)) ;
imageSQ(tt, tfrticC*basicTF.fs, rtfrC, 0.995); axis xy; colormap(1-gray);
hold on ; plot(t, ones(size(t)), 'r') ;
plot(t(1:basicTF.hop:end), IF2(1:20*basicTF.hop:end), 'r') ;
xlabel('time (s)'); ylabel('frequency (Hz)');

pause



%===================================
	% reconstruct component
alpha = Hz * (tfrsqtic(2)-tfrsqtic(1)) ;
RR = floor(0.8./alpha);
xhat = zeros(size(x));
x2hat = zeros(size(x));
	
		%% (cheat) extract curves
tmp_tfrsq1 = zeros(31, size(tfrsq,2));
idx1 = floor((clean{1}.IF - tfrsqtic(1))./alpha)+1;
tmp_tfrsq2 = zeros(31, size(tfrsq,2));
idx2 = floor((clean{2}.IF - tfrsqtic(1))./alpha)+1;

for ii = 1:size(tfrsq,2)
    tmp_tfrsq1(:, ii) = abs(tfrsq(idx1(ii) - 15:idx1(ii) + 15, ii));
    tmp_tfrsq2(:, ii) = abs(tfrsq(idx2(ii) - 15:idx2(ii) + 15, ii));
end

c1 = CurveExt_M(tmp_tfrsq1', 10);
c1 = max(1, idx1-15)-1+c1;
c2 = CurveExt_M(tmp_tfrsq2', 10);
c2 = max(1, idx2-15)-1+c2;


for oo = 1: length(x)
    xhat(oo) = h((opts.WinLen+1)/2)*sum(tfrsq(max(1,c1(oo)-RR): min(length(tfrsqtic),c1(oo)+RR), oo))*alpha;
    x2hat(oo) = h((opts.WinLen+1)/2)*sum(tfrsq(max(1,c2(oo)-RR): min(length(tfrsqtic),c2(oo)+RR), oo))*alpha;
end



%=================
	%% estimate shape

ifhat = tfrsqtic(c1);
if2hat = tfrsqtic(c2);
phihat0 = cumsum(2*pi*ifhat) ./ Hz ;
phi2hat0 = cumsum(2*pi*if2hat) ./ Hz ;

%phihat0 = 2*pi*data.clean.phase(1,:);
%phi2hat0 = 2*pi*data.clean.phase(2,:);

ahat = abs(xhat);
a2hat = abs(x2hat);
phihat = phase(xhat ./ ahat);
phi2hat = phase(x2hat ./ a2hat);


idx = find(isnan(phihat)) ;
if length(idx)
    idx = setdiff([1:length(c1)], idx);
    phihat = interp1(idx,phihat(idx),[1:1024]');
end

idx = find(isnan(phi2hat))
if length(idx)
    idx = setdiff([1:length(c1)], idx);
    phi2hat = interp1(idx,phi2hat(idx),[1:1024]');
end


%===============
    %% better estiamte AM
%[A] = estAM(xm, fm)

lambda = 1e3 ;
C = [diag(ones(n,1)) diag(ones(n,1))] ;
Phi = diag([cos(phihat); cos(phi2hat)]) ;

D0 = -diag(ones(n, 1));

for ii = 1 : n-2
    D0(ii, ii+1) = 1;
end

D0(n-1, n) = 1;



%{
D0 = zeros(n) ;
for jj = 2:n-1
    D0(jj, jj-1) = 1 ; D0(jj,jj) = -2 ; D0(jj, jj+1) = 1;
end
D0(n, 1) = 1 ; D0(n, n) = -2 ; D0(n, n-1) = 1 ;
D0(1, 1) = -2 ; D0(1, 2) = 1 ; D0(1, n) = 1 ;
%}

D = zeros(2*n) ;
D(1:n,1:n) = D0 ; D(n+1:2*n, n+1:2*n) = D0 ;
A = xm' * C * Phi * inv(Phi' * C' * C * Phi + lambda * D'*D) ;
ahat = A(1:n)' ;
a2hat = A(n+1:2*n)' ;

%===============
ahat = ahat(50:end-49);
a2hat = a2hat(50:end-49);
phihat = phihat(50:end-49);
phi2hat = phi2hat(50:end-49);


    %% extract the shape
    %% suppose there are 5 components in the shape
K = 5;

    %% the shape estimator
c = zeros(K*2*2, length(phihat));
for ii = 1:K
    c((ii-1)*2+1,:) = ahat .* cos(ii*phihat);
    c((ii-1)*2+2,:) = ahat .* sin(ii*phihat);
    c(K*2+(ii-1)*2+1,:) = a2hat .* cos(ii*phi2hat);
    c(K*2+(ii-1)*2+2,:) = a2hat .* sin(ii*phi2hat);
end

d = zeros(K*2*2, length(phihat));
for ii = 1:K
    d((ii-1)*2+1,:) = cos(ii*phihat);
    d((ii-1)*2+2,:) = sin(ii*phihat);
    d(K*2+(ii-1)*2+1,:) = cos(ii*phi2hat);
    d(K*2+(ii-1)*2+2,:) = sin(ii*phi2hat);
end

hatshat = (xm(50:end-49)' * c') * inv(c*c');

s1hat = hatshat(1:K*2) * c(1:K*2, :);
s2hat = hatshat(K*2+1:end) * c(K*2+1:end, :);

MAX = max(xm(50:end-49))*0.7;
MIN = min(xm(50:end-49))*0.7;

subplot(2,1,1); hold off;
plot(t, clean{1}.component); hold on; plot(t(50:end-49), s1hat,'r'); axis tight; set(gca,'fontsize',24);

subplot(2,1,2); hold off;
plot(t, clean{2}.component); hold on; plot(t(50:end-49), s2hat,'r'); axis tight; set(gca,'fontsize',24);

%===============


