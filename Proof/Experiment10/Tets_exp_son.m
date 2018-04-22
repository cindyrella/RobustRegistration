clear all;
n =100000;
X1 = randn(n,2);
X2 = randn(n,2);
rho1 = (X1.^2)*ones(2,1);
rho2 = (X2.^2)*ones(2,1);
theta1 = pi/3;
theta2 = pi/6;

den = 2*sqrt((sin(theta1/2))^2*rho1 + (sin(theta2/2))^2*rho2);
exp1 = sum(rho1./den)/n;
exp2 = sum(rho2./den)/n;

ideal1 = sqrt(pi/2)/2*(1/(sin(theta1/2)+sin(theta2/2))^2)*(sin(theta1/2)+2*sin(theta2/2));
ideal2 = sqrt(pi/2)/2*(1/(sin(theta1/2)+sin(theta2/2))^2)*(sin(theta2/2)+2*sin(theta1/2));