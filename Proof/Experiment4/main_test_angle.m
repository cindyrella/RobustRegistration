clear all;
%% Generate data
p = 0.5;
psi = pi/2;
n = 1000;
X = randn(n,3);
P = binornd(1,p,n,1);
Y = randn(n,3).*P+X.*(1-P);

XY = diag(X*Y');
YcX = cross(Y,X);
XmY = diag((X-Y)*(X-Y)');

%% Generate angles
m = 20;
thetavec = linspace(0,2*pi,2*m);
phivec = linspace(0,pi,m);
[ph, th] = meshgrid(phivec,thetavec);

R = ones(size(th)); 

b = R.*sin(ph).*cos(th);
c = R.*sin(ph).*sin(th);
d = R.*cos(ph);

Col1 = zeros(numel(thetavec),numel(phivec));
figure()
hold on;
for id_t = 1:numel(thetavec)
    id_t
    for id_p = 1:numel(phivec)
        w = [b(id_t,id_p),c(id_t,id_p),d(id_t,id_p)]';
        test_angle(w,X,Y,XY,YcX,XmY);
    end
end
plot(thetavec,thetavec);
hold off;



