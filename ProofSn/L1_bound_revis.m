m = 100;
theta = linspace(0,pi/2,m);
for n = 1:10
    for l=[0.1,0.2,0.3]
        s = sqrt(pi)./4./(1+(sqrt(n)-1)*l./sin(theta));
        s = 2*asin(s);
        s = (theta<s).*s;
        plot(theta,s);
        hold on
    end
end
        