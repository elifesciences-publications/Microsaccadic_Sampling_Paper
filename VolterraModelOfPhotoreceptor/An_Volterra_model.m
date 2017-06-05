% simple Volterra series system identification (only zero & first order)
% May 2010
% Written by An Dau

function [k0, k1]  = Vol_model(x, y, kernellength, windowlength)
% kernel and window length are both in data points
% x and y are input and output of the training data
R = kernellength;
N = windowlength;

% matrix Y
Y(N,1) = 0;
for k = 1:N
    Y(k,1) = y(N+R-k);
end
% matrix X1
X1(N,(R+1))=0;
for i = 1:N
    X1(i,1) = 1;
    for j = 2:(R+1)
        X1(i,j) = x(N+R+2-i-j);
    end
end

K1 = pinv(X1)*Y;

% output values
k0 = K1(1);
k1 = K1(2:(R+1));
plot(k1)