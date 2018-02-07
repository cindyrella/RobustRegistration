d=3;
n=4;

A = randn(d*n);
B = triu(A)+triu(A)';

for i = 1:n
    B(((i-1)*d+1):(i*d),((i-1)*d+1):(i*d))=eye(d);
end

eig(B)