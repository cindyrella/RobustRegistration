%% Plot cosine part
m = 100;
theta = linspace(0,2*pi,m);
figure()
plot(theta,(1-cos(theta)));
hold on;
plot(theta,sqrt(2*(1-cos(theta))));
hold off;

figure()
plot(theta,-sin(theta));
hold on;
plot(theta,1./sqrt(2*(1-cos(theta))).*sin(theta));
hold off;
