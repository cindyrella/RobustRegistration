%n = 1000000;
n = 1000
thetaL1 = zeros(1,n);
thetaL2 = zeros(1,n);
% for j = 1:10
% p = 0.1*j;
p = 0.5;
for i = 1:n
m_tot = 100;
m = ceil(p*m_tot);
m_clean = m_tot - m;
X = randn(m,2);
Y = randn(m,2);
X_clean = randn(m_clean,2);
A = [0,1;-1,0];
theta_2 = atan(trace(A*X'*Y)/(trace(X'*Y)+trace(X_clean'*X_clean)));
if (trace(X'*Y)+trace(X_clean'*X_clean)) < 0
    theta_2 = theta_2 + pi;
end
thetaL2(i) = theta_2; 
thetaL1(i) = approx_angle([X;X_clean],[Y;X_clean]);
end

figure()
hist(thetaL2,30);
hold on;
hist(thetaL1,30);
hold off
%end