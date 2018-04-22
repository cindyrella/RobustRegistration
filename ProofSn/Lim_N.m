N = 1:10;
alpha = [0.15,0.1,0.05,0.01];
min_N_L1 =@(n,a) (64/pi) .*(log(n)-log(a)).*(n+sqrt((pi/2).*n)).^2;
min_N_L2 =@(n,a) 8 .*(log(n)-log(a)).*(n+1).^2;

colors = jet(numel(alpha));
for a = 1:numel(alpha)
    plot(N,min_N_L1(N,alpha(a)),'-','Color',colors(a,:));
    hold on;
    plot(N,min_N_L2(N,alpha(a)),'--','Color',colors(a,:));
end