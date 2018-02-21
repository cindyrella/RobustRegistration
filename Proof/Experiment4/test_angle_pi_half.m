%% Analize unitary rotation
clear all;
%% Generate data
p = 0.8;
psi = pi/2;
n = 1000;
X = randn(n,3);
P = binornd(1,p,n,1);
Y = randn(n,3).*P+X.*(1-P);

XY = (X'*Y+Y'*X)/2;
YcX = sum(cross(Y,X));
[~,S,V] = svd(X,0);
U = [V(:,3),V(:,1),V(:,2)];
%% Generate angles
m = 100;
thetavec = linspace(0,2*pi,2*m);
phivec = linspace(0,pi,m);
[ph, th] = meshgrid(phivec,thetavec);

R = ones(size(th)); 

b = R.*sin(ph).*cos(th);
c = R.*sin(ph).*sin(th);
d = R.*cos(ph);

Col1 = zeros(numel(thetavec),numel(phivec));
for id_t = 1:numel(thetavec)
    id_t
    for id_p = 1:numel(phivec)
        w = U*([b(id_t,id_p),c(id_t,id_p),d(id_t,id_p)]');
        
        Col1(id_t,id_p) = atan(YcX*w/(trace(XY)-w'*(XY)*w))/pi*180;
    end
end
%% Plot Norms
s = surf(b,c,d,Col1);
hold on;
s.EdgeColor = 'none';
scatter3([1,0,0],[0,1,0],[0,0,1],'*k');
%quiver3(b,c,d,Grad1(:,:,1),Grad1(:,:,2),Grad1(:,:,3))
hold off;

cvx_begin sdp
    variable w(3) 
    variable theta(1)
    minimize -1*sin(theta)*w'*YcX+(1-cos(theta))*(trace(XY)-w'XY*w)
    w'*w ==1
cvx_end
