%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  *** Examples demonstrating the use of the code ***
%  
%
%   (C) Martin Slawski, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%***[1] example in 3D, cf. right plot of Fig. 1 in the NIPS 2013 paper ***%

T= [0 0 0; 1 1 0; 0 1 1]';

A = [1/3 1/3 1/3; 
        3/4 1/8 1/8; 
        1/6 2/3 1/6;
        3/20 7/10 3/20;
        1/12 1/12 5/6;
        1/20  4/5  1/20;
        1/4    1/2  1/4; 
        3/10  2/5   3/10;
        2/5    3/10 3/10;
        3/10  3/10 2/5]';
    
D = T * A;

opt0 = opt_Integerfac_findvert();

[That, Ahat, status] = Integerfac_findvert_cpp(D, 3, [0 1], opt0);

%***[2]  random example ***% 

T = double(rand(100, 10) > 0.5);

A = randn(10, 50);

A(end,:) = 1 - sum(A(1:end-1,:));
% all(sum(A) == 1)

D = T * A;

% we have to remove the non-negativity constraints on A,
% which is the default
opt_anysign = opt_Integerfac_findvert('nonnegative', false);

[That, Ahat, status] = Integerfac_findvert_cpp(D, 10, [0 1], opt_anysign);

%max(min(dist_euclidean(That', T')))
%max(min(dist_euclidean(Ahat, A)))

% we can also drop the sum-to-one constraint on A
% by modifying the options:

opt_linear = opt_Integerfac_findvert('nonnegative', false, 'affine', false);

[That, Ahat, status] = Integerfac_findvert_cpp(D, 10, [0 1], opt_linear);

% we can deal with {-1,1}-matrices as well:

T = 2*T - 1; % {-1,1}
D = T * A;

[That, Ahat, status] = Integerfac_findvert_cpp(D, 10, [-1 1], opt_anysign);

%%% ***[3] presence of noise *** %%%

E = 0.1 * randn(size(D));
D = D + E;

[That, Ahat, status] = Integerfac_findvert_cpp(D, 10, [-1 1], opt_anysign);

%max(min(dist_euclidean(That', T')))
%max(min(dist_euclidean(Ahat, A)))

% heavier noise

D = That * Ahat + 10*E;

[That, Ahat, status] = Integerfac_findvert_cpp(D, 10, [-1 1], opt_anysign);

% error can often be reduced by aggregating over multiple
% selections of rows

opt_aggr =  opt_anysign;
opt_aggr.aggr = 'bestsingle';

[That, Ahat, status] = Integerfac_findvert_cpp(D, 10, [-1 1], opt_aggr);

%max(min(dist_euclidean(That', T')))
%max(min(dist_euclidean(Ahat, A)))

%%% *** [4] Ternary matrix factorization *** %%%

% takes values in {0, 0.5, 1}
T = projectV(rand(100, 10), [0 0.5 1]);
D = T * A;

[That, Ahat, status] = Integerfac_findvert_cpp(D, 10, [0 0.5 1], opt_anysign);

%max(min(dist_euclidean(That', T')))
%max(min(dist_euclidean(Ahat, A)))

%%% *** End of file *** %%%%%%%%%%%%%%%%%%
