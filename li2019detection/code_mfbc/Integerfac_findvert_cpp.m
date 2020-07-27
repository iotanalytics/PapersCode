function [T, A, status] = Integerfac_findvert_cpp(D, r, V, opt)
%*************************************************************************% 
%
% -- Input --
%
%   D      - input data matrix
%   V      - finite set of values the components are allowed to take
%   r       - number of components
%   opt    - various options, see Integerfac_findvert_opt
%            Default options are used in case that this argument is omitted.
%
%
% -- Output --
%
%   T       -  matrix of components
%
%   A       -  coefficient matrix
%
%   status  - a struct containing ()additional status information
%
%
%   Reference:
%
%   Martin Slawski, Matthias Hein and Pavlo Lutsik
%   'Matrix Factorization with Binary Components' 
%    NIPS 2013
%
%    This is an implementation of an extended version of Algorithm 3
%    in the above paper. 
%
%
%   (C)  Martin Slawski, 2014  
%
%*************************************************************************%

if nargin < 4
    opt = opt_Integerfac_findvert(); % initialize with default options
end

if nargin < 3 || isempty(V)
    V = [0 1];
end

if ~issorted(V)
    error('V must be sorted in ascending order')
end

if ismember(0, V) &&  ~opt.affine % in order to throw out the all-zero-vector
    zeroprob = true;
else
    zeroprob = false;
end

[m, n] = size(D);

if  opt.affine 
    if ~(r-1 < min(m, n))
        error('r-1 must be smaller than min(m, n)')
    end
else
    if ~(r  < min(m, n))
        error('r must be smaller than min(m, n)')
    end
end

cardV = length(V);
myeps = 1E-10;

%

if opt.affine
    meanD = mean(D,2);
    E = D - meanD*ones(1,n);
    k = r-1;
else
    meanD = zeros(m, 1); % to simplify the code
    E = D;
    k = r;
end

%
[U, sigma, ~] = svds(E, k);
if sigma(end) < sqrt(myeps)
    error('Data matrix does not match the specified number of components r: one of the top singular values is (close) to zero')
end
ixx = 1:k;
TP=[]; % TP is the matrix for the found topics (vertices)

All = 1:m;

if strcmp(opt.aggr, 'no')
    nsamples_basic = 1; nsamples_extra = 0;
else
    nreplace = round(opt.replace * k);
    nreplace = max(1, nreplace); % at least one coordinate has to be replaced
    nsamples_basic = max(floor((m - k)/nreplace), 1);
    
    if opt.nsamples == 0
        nsamples_extra = 0;
    else
        if opt.nsamples > 0
            nsamples_tot = opt.nsamples;
            if nsamples_tot > nsamples_basic
                nsamples_extra = nsamples_tot - nsamples_basic;
            else
                nsamples_basic = min(opt.nsamples, nsamples_basic); nsamples_extra = 0;
            end
        else
            nsamples_extra = -(opt.nsamples);
        end
    end
end
nsamples = nsamples_basic + nsamples_extra;

%%%

if  opt.nonnegative
    if opt.affine
        A0 = randsplxmat(r,n);
    else
        A0 = rand(r,n);
    end
end

%%%

if(k>2)
    if cardV > 2
        rest=max(k-opt.chunksize,1);
        PossOver_0 = transpose(combn(V, rest));
        Poss_0 = transpose(combn(V, k-rest));
    end
else
    Poss_0 = transpose(combn(V, k));
end

%%%

for i=1:nsamples
    if  i <= nsamples_basic
        [~,kxx] = licols(U(All,ixx)');
        jxx=All(kxx);
        if ~strcmp(opt.aggr, 'no'); All(kxx(1:nreplace))=[]; end
        
        if opt.verbose; disp(['Round: ',num2str(i),' - Selecting: ',num2str(jxx)]); end
    else % choose extra coordinates by random selection
        rpm = randperm(m);
        jxx = rpm(1:k);
        if opt.verbose; disp(['Round: ',num2str(i),' - Selecting: ',num2str(jxx)]); end
    end
    
    FD=U(jxx,ixx);
    v = FD \ eye(k);
    TT = U(:,ixx)*v;
    Indices = 1:m;  % we only have to test for the rows (Indices) which were not selected above (jxx)
    Indices(jxx) = [];
    
    meanSmall = TT * meanD(jxx);
    
    if(k>2)
        
        if cardV == 2
            
            scale = V(2) - V(1);
            shift = V(1);
            
            TT_shift = shift * sum(TT, 2);
            
            [TPP, TPPPattern, dist] = getVertices(scale * TT(Indices,:)', TT_shift(Indices) + meanD(Indices)-meanSmall(Indices), V);
            
            if zeroprob
                
                % sort out the zero vector
                TPP = TPP(:,2:end);
                TPPPattern = TPPPattern(:,2:end)  *  scale + shift;
            
             else
                 
                 [~, ix_so_dist] = sort(dist, 'ascend');
                TPP = TPP(:,ix_so_dist(1:r));
                TPPPattern = TPPPattern(:, ix_so_dist(1:r))  *  scale + shift;
                    
            end
            
           
            
            TP(Indices,  ((i-1)*r+1) : i*r) = TPP;
            TP(jxx,  ((i-1)*r+1) : i*r) = TPPPattern * scale + shift;
            
        else
            
            PossOver= PossOver_0 - repmat(meanD(jxx(1:rest)),1,cardV^rest);
            Poss=  Poss_0 - repmat(meanD(jxx(rest+1:k)),1,cardV^(k-rest));
            
            num = size(PossOver,1); nump=num+1;
            Poss = [Poss;ones(1,size(Poss,2))]; % we append a one to the vector to integrate the fixed part
            
            TT1=TT(Indices,1:num); TT2=TT(Indices,nump:end);   % thus we also select just the rows Indices from T (the rest is the identity matrix) and
            % divide it up into the fixed component (rest) and the chunk
            
            for kk=1:size(PossOver,2)
                
                vec=PossOver(:,kk);
                TPP1 = TT1*vec; % these are the first components
                TPP = [TT2,TPP1]*Poss+repmat(meanD(Indices),1,size(Poss,2));   % these are the first offset rows of the potential topics
                TTT=abs(TPP - projectV(TPP, V));
                if zeroprob && kk==1
                    TTT = TTT(:,2:end); % get rid of the all-zero-vector
                end
                [sdist,icc]=sort(sum(TTT.^2),'ascend');
                icc = icc + (zeroprob && kk==1);
                
                if(kk==1)
                    X1=TPP(:,icc(1:r));
                    X2 = [repmat(vec+meanD(jxx(1:rest)),1,r);Poss(1:end-1,icc(1:r))+repmat(meanD(jxx(rest+1:k)),1,r)];
                    if(i==1)
                        TP=zeros(size(D,1),size(X1,2));
                    end
                    TP(Indices,((i-1)*r+1) : i*r)=X1;
                    TP(jxx,((i-1)*r+1): i*r)=X2;
                    dT(1:r)=sdist(1:r);
                    tjx=icc(1:r);
                else
                    mixdT=[dT,sdist(1:r)]; % join new and old distances
                    [smixdT,idd]=sort(mixdT,'ascend');
                    thresh=smixdT(r); % this is the new maximal distance
                    
                    freeslots=find(dT> thresh);
                    tjx = icc(1:length(freeslots));
                    for ii=1:length(freeslots)
                        X1=TPP(:,tjx(ii));
                        X2=[vec+meanD(jxx(1:rest));Poss(1:end-1,tjx(ii))+meanD(jxx(rest+1:k))];
                        TP(Indices,(i-1)*r+freeslots(ii))=X1;
                        TP(jxx,(i-1)*r+freeslots(ii))=X2;
                    end
                    dT(freeslots)=sdist(1:length(freeslots));
                end
            end
            
        end
    else
        Poss = Poss_0 - repmat(meanD(jxx(1:k)),1,cardV^k);
        TPP = TT*Poss + repmat(meanD, 1, cardV^k);
        TTT=abs(TPP - projectV(TPP, V));
        [~,icc]=sort(sum(TTT.^2),'ascend');
        TP(:,((i-1)*r+1) : i*r)=TPP(:,icc((1+zeroprob):(r+zeroprob)));
    end
    
    if opt.verbose; disp(['Single Fit - Iteration: ',num2str(i)]); end
    TP = projectV(TP, V);
    T2=TP(:,((i-1)*r+1):i*r);
    if opt.nonnegative
        G = T2' * T2; W = T2' * D;
        if opt.affine
            [A2, ~, ~] = mexQuadSimplex(G,W,A0, myeps);
        else
            [A2, ~, ~] = mexQuadNN(G,W,A0, myeps);
        end
    else
        if opt.affine
            A2 = affls(T2, D);
        else
            A2 = T2 \ D; % least squares
        end
    end
    singleerr(i)= sum(sum(abs(D-T2*A2).^2));
end

[bestsingleerr,mix]=sort(singleerr,'ascend');
bestsingleerr=bestsingleerr(1);
if opt.verbose; disp(['Best single error: ',num2str(bestsingleerr)]); end

if ismember(opt.aggr, {'no', 'bestsingle'})
    T = TP(:,((mix(1)-1)*r+1):mix(1)*r);
    if opt.nonnegative
        G = T' * T; W = T' * D;
        if opt.affine
            [A, ~, ~] = mexQuadSimplex(G,W,A0, myeps);
        else
            [A, ~, ~] = mexQuadNN(G,W,A0, myeps);
        end
    else
        if opt.affine
            A = affls(T, D);
        else
            A = T \ D;
        end
    end
    status.err = sum(sum(abs(D-T*A).^2)); status.nsamples = nsamples;
    return;
end

% merge the best candidate sets
for i=1:min(nsamples, opt.naggsets)
    TP2(:,((i-1)*r+1):i*r)=TP(:,((mix(i)-1)*r+1):mix(i)*r);
end
TP=TP2;

TC(:,1)=TP(:,1);
counter=2;
for i=2:size(TP,2)
    if(dist_euclidean(TC',TP(:,i)')>1E-3*m)
        TC(:,counter)=TP(:,i);
        counter=counter+1;
    end
end
clear TP;
if opt.verbose; disp(['Number of topics after duplicate checking: ',num2str(size(TC,2))]); end

%%%%

if opt.nonnegative
    if opt.affine
        T2=TC; G = T2' * T2; A0 = randsplxmat(size(T2,2),n); W = T2' * D;
        [A2, ~, ~] = mexQuadSimplex(G,W,A0, myeps);
        err=sum(sum(abs(D-T2*A2).^2));
        if opt.verbose; disp(['With ',num2str(size(TC,2)),' topics - err: ',num2str(err)]); end
        
        while(size(TC,2)>r)
            clear err;
            for i=1:size(TC,2)
                itt=setdiff(1:size(TC,2),i);
                T2=TC(:,itt);
                if(rem(i,10)==0), if opt.verbose; disp(['To do: ',num2str(size(TC,2)),' - done: ',num2str(i)]); end; end
                G = T2' * T2; A0 = randsplxmat(size(T2,2),n); W = T2' * D;
                [A2, ~, ~] = mexQuadSimplex(G,W,A0, myeps);
                err(i)=sum(sum(abs(D-T2*A2).^2));
            end
            if(size(TC,2)<4*r)
                [~,disc]=min(err);
                if opt.verbose; disp(['Current size: ',num2str(size(TC,2)),' - Minimal error: ',num2str(min(err)),' - discarding: ',num2str(disc)]); end
                itt=setdiff(1:size(TC,2),disc);
                TC=TC(:,itt);
            else
                dec=5;
                if(size(TC,2)>k+50), dec=20; end
                if(size(TC,2)>k+100), dec=50; end
                if(size(TC,2)>k+200), dec=100; end
                if(size(TC,2)>k+2000), dec=1000; end
                [serr,iuu]=sort(err,'ascend');
                if opt.verbose; disp(['Minimal error: ',num2str(serr(1)),' Error ',num2str(dec),': ',num2str(serr(dec)),' - discarding: ',num2str(iuu(1:dec))]); end
                itt=setdiff(1:size(TC,2),iuu(1:dec));
                TC=TC(:,itt);
            end
        end
        
        % final fit
        T2=TC; G = T2' * T2; A0 = randsplxmat(size(T2,2),n); W = T2' * D;
        [A2, ~, ~] = mexQuadSimplex(G,W,A0, myeps);
        
    else
        
        T2=TC; G = T2' * T2; A0 = rand(size(T2,2),n); W = T2' * D;
        [A2, ~, ~] = mexQuadNN(G,W,A0, myeps);
        err=sum(sum(abs(D-T2*A2).^2));
        if opt.verbose; disp(['With ',num2str(size(TC,2)),' topics - err: ',num2str(err)]); end
        
        while(size(TC,2)>r)
            clear err;
            for i=1:size(TC,2)
                itt=setdiff(1:size(TC,2),i);
                T2=TC(:,itt);
                if(rem(i,10)==0), if opt.verbose; disp(['To do: ',num2str(size(TC,2)),' - done: ',num2str(i)]); end; end
                G = T2' * T2; A0 = rand(size(T2,2),n); W = T2' * D;
                [A2, ~, ~] = mexQuadNN(G,W,A0, myeps);
                err(i)=sum(sum(abs(D-T2*A2).^2));
            end
            if(size(TC,2)<4*r)
                [~,disc]=min(err);
                if opt.verbose; disp(['Current size: ',num2str(size(TC,2)),' - Minimal error: ',num2str(min(err)),' - discarding: ',num2str(disc)]); end
                itt=setdiff(1:size(TC,2),disc);
                TC=TC(:,itt);
            else
                dec=5;
                if(size(TC,2)>k+50), dec=20; end
                if(size(TC,2)>k+100), dec=50; end
                if(size(TC,2)>k+200), dec=100; end
                if(size(TC,2)>k+2000), dec=1000; end
                [serr,iuu]=sort(err,'ascend');
                if opt.verbose; disp(['Minimal error: ',num2str(serr(1)),' Error ',num2str(dec),': ',num2str(serr(dec)),' - discarding: ',num2str(iuu(1:dec))]); end
                itt=setdiff(1:size(TC,2),iuu(1:dec));
                TC=TC(:,itt);
            end
        end
        
        % final fit
        T2=TC; G = T2' * T2; A0 = rand(size(T2,2),n); W = T2' * D;
        [A2, ~, ~] = mexQuadNN(G,W,A0, myeps);
        
    end
    
else
    
    if opt.affine
        T2=TC;
        A2 = affls(T2, D);
        err=sum(sum(abs(D-T2*A2).^2));
        if opt.verbose; disp(['With ',num2str(size(TC,2)),' topics - err: ',num2str(err)]); end
        
        while(size(TC,2)>r)
            clear err;
            for i=1:size(TC,2)
                itt=setdiff(1:size(TC,2),i);
                T2=TC(:,itt);
                if(rem(i,10)==0), if opt.verbose; disp(['To do: ',num2str(size(TC,2)),' - done: ',num2str(i)]); end; end
                A2 = affls(T2, D);
                err(i)=sum(sum(abs(D-T2*A2).^2));
            end
            if(size(TC,2)<4*r)
                [~,disc]=min(err);
                if opt.verbose; disp(['Current size: ',num2str(size(TC,2)),' - Minimal error: ',num2str(min(err)),' - discarding: ',num2str(disc)]); end
                itt=setdiff(1:size(TC,2),disc);
                TC=TC(:,itt);
            else
                dec=5;
                if(size(TC,2)>k+50), dec=20; end
                if(size(TC,2)>k+100), dec=50; end
                if(size(TC,2)>k+200), dec=100; end
                if(size(TC,2)>k+2000), dec=1000; end
                [serr,iuu]=sort(err,'ascend');
                if opt.verbose; disp(['Minimal error: ',num2str(serr(1)),' Error ',num2str(dec),': ',num2str(serr(dec)),' - discarding: ',num2str(iuu(1:dec))]); end
                itt=setdiff(1:size(TC,2),iuu(1:dec));
                TC=TC(:,itt);
            end
        end
        
        % final fit
        T2 = TC;
        A2 = affls(T2, D);
        
    else
        
        T2 = TC;
        A2 = T2 \ D;
        err=sum(sum(abs(D-T2*A2).^2));
        if opt.verbose; disp(['With ',num2str(size(TC,2)),' topics - err: ',num2str(err)]); end
        
        while(size(TC,2)>r)
            clear err;
            for i=1:size(TC,2)
                itt=setdiff(1:size(TC,2),i);
                T2=TC(:,itt);
                if(rem(i,10)==0), if opt.verbose; disp(['To do: ',num2str(size(TC,2)),' - done: ',num2str(i)]); end; end
                A2 = T2 \ D;
                err(i)=sum(sum(abs(D-T2*A2).^2));
            end
            if(size(TC,2)<4*r)
                [~,disc]=min(err);
                if opt.verbose; disp(['Current size: ',num2str(size(TC,2)),' - Minimal error: ',num2str(min(err)),' - discarding: ',num2str(disc)]); end
                itt=setdiff(1:size(TC,2),disc);
                TC=TC(:,itt);
            else
                dec=5;
                if(size(TC,2)>k+50), dec=20; end
                if(size(TC,2)>k+100), dec=50; end
                if(size(TC,2)>k+200), dec=100; end
                if(size(TC,2)>k+2000), dec=1000; end
                [serr,iuu]=sort(err,'ascend');
                if opt.verbose; disp(['Minimal error: ',num2str(serr(1)),' Error ',num2str(dec),': ',num2str(serr(dec)),' - discarding: ',num2str(iuu(1:dec))]); end
                itt=setdiff(1:size(TC,2),iuu(1:dec));
                TC=TC(:,itt);
            end
        end
        
        % final fit
        T2 = TC;
        A2 = T2 \ D;
        
    end
    
end

err=sum(sum(abs(D-T2*A2).^2));
T=T2; A=A2;
if opt.verbose; disp(['With ',num2str(size(TC,2)),' topics - err: ',num2str(err)]); end
status.nsamples = nsamples;
status.err = err;

return;