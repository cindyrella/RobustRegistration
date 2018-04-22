%% Test that the variance is 1/32 of the L1 noise case
n = 1000000;
Y = randn(4,n);
S1 = 0;
S2 = 0;
for i = 1:n
    X1 = Y(1:2,i);
    X2 = Y(3:4,i);
    cross = X1(1)*X2(2) - X1(2)*X2(1);
    dot = X1'* X2;
    total = cross/sqrt(X1'*X1 + X2'*X2 - 2*dot);
    S1 = S1+total^2;
    S2 = S2+total;
end
S1/(n-1)
S2/n