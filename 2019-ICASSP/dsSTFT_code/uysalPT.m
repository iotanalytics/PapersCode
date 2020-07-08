function [Q_loc,Q_value,R_loc,R_value,S_loc,S_value] = uysalPT(ecg,fs)
x1 = ecg;
N = length (x1);       % Signal length
t = [0:N-1]/fs;        % time index

%%%Cancellation DC drift & normalization
x1 = x1 - mean (x1 );    % cancel DC conponents
x1 = x1/ max( abs(x1 )); % normalize to one

%%%Low pass filtering
% LPF (1-z^-6)^2/(1-z^-1)^2
b=[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a=[1 -2 1];
 
 
h_LP=filter(b,a,[1 zeros(1,12)]); % transfer function of LPF
 
x2 = conv (x1 ,h_LP);
%x2 = x2 (6+[1: N]); %cancle delay
x2 = x2/ max( abs(x2 )); % normalize , for convenience .

%%%High pass filtering
% HPF = Allpass-(Lowpass) = z^-16-[(1-z^-32)/(1-z^-1)]
b = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a = [1 -1];
 
h_HP=filter(b,a,[1 zeros(1,32)]); % impulse response iof HPF
 
x3 = conv (x2 ,h_HP);
%x3 = x3 (16+[1: N]); %cancle delay
x3 = x3/ max( abs(x3 ));

%%%Derivative filter
% Make impulse response
h = [-1 -2 0 2 1]/8;
% Apply filter
x4 = conv (x3 ,h);
x4 = x4 (2+[1: N]);
x4 = x4/ max( abs(x4 ));

%%%Squaring
x5 = x4 .^2;
x5 = x5/ max( abs(x5 ));

%%%Moving window integration
% Make impulse response
h = ones (1 ,31)/31;
Delay = 15; % Delay in samples
 
% Apply filter
x6 = conv (x5 ,h);
x6 = x6 (15+[1: N]);
x6 = x6/ max( abs(x6 ));

%%%Find QRS points
max_h = max(x6);
thresh = mean (x6 );
poss_reg =(x6>thresh*max_h)';

left = find(diff([0 poss_reg])==1);
right = find(diff([poss_reg 0])==-1);
 
left=left-(6+16);  % cancle delay because of LP and HP
right=right-(6+16);% cancle delay because of LP and HP

for j = 1:length(left)
    if left(j) < 0; 
        left(j) = 1;
    end
end

for k = 1:length(right)
    if right(k) < 0; 
        right(k) = 1;
    end
end

for i=1:length(left)
    [R_value(i) R_loc(i)] = max( x1(left(i):right(i)) );
    R_loc(i) = R_loc(i)-1+left(i); % add offset
 
    [Q_value(i) Q_loc(i)] = min( x1(left(i):R_loc(i)) );
    Q_loc(i) = Q_loc(i)-1+left(i); % add offset
 
    [S_value(i) S_loc(i)] = min( x1(left(i):right(i)) );
    S_loc(i) = S_loc(i)-1+left(i); % add offset
 
end
 
% there is no selective wave
Q_loc=Q_loc(find(Q_loc~=0));
R_loc=R_loc(find(R_loc~=0));
S_loc=S_loc(find(S_loc~=0));
