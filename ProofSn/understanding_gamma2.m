ntheta = 100
n=2


theta = linspace(0,pi,100);
delta = linspace(0,pi,10);
for i = 1:numel(delta)
    theta2 = min(theta+delta(i),pi);
    g1 = max(sqrt(1-cos(theta)),sqrt(1-cos(theta2)));
    g2 = max(sqrt(1+cos(theta)),sqrt(1+cos(theta2)));
    g3 = max(max(abs(sin(theta)),abs(sin(theta2))),abs(sin((theta+theta2)/2)));
    g3p =max(abs(sin(theta)),abs(sin(theta2)));
    
%     plot(theta,g1.*g2./g3);
%     hold on;
%     plot(theta,g1.*g2./g3p,'--');
% plot(theta,g1.*g2+2/n*g3)
% hold on;
% plot(theta,g1.*g2+2/n*g3p)
plot(theta,min(1,g1.*g2+2/n*g3p).^2/max(2,g1.*g2+2/n*g3p)./g3p)
end