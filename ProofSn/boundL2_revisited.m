n = 1:10;
p = 0:0.05:0.5;
A = [0.15,0.1,0.05,0.1];
[N,P]=meshgrid(n,p);
for alpha = A
    casa = (log(N)-log(alpha)).*2.*N.^2.*(P+4*(1-P).*(1+1./N).^2)./(1-P).^2;
    casa = sqrt(casa);
    figure()
    imagesc(casa);
end