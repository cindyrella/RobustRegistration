%% Analize unitary rotation

%% Generate data
n = 1000;
X = randn(n,3);
[~,S,V] = svd(X,0);
%% Generate angles
m = 20;
thetavec = linspace(0,2*pi,2*m);
phivec = linspace(0,pi,m);
[ph, th] = meshgrid(phivec,thetavec);

R = ones(size(th)); 

b = R.*sin(ph).*cos(th);
c = R.*sin(ph).*sin(th);
d = R.*cos(ph);

%% Add differences

Col1 = zeros(numel(thetavec),numel(phivec));
Col2 = zeros(numel(thetavec),numel(phivec));
for id_t = 1:numel(thetavec)
    for id_p = 1:numel(phivec)
        accum1 = 0;
        accum2 = 0;
        w = [b(id_t,id_p),c(id_t,id_p),d(id_t,id_p)]';
        for id_X = 1:n
            x = X(id_X,:)';
            accum1 = accum1 + sqrt(x'*x - (x'*w)^2);
            accum2 = accum2 + x'*x - (x'*w)^2;
        end    
        Col1(id_t,id_p) = accum1;
        Col2(id_t,id_p) = accum2;
    end
end
s = surf(b,c,d,Col1);
hold on;
s.EdgeColor = 'none';
scatter3(V(1,:),V(2,:),V(3,:),'*k');
hold off;
figure()
s = surf(b,c,d,Col2);
hold on;
s.EdgeColor = 'none';
scatter3(V(1,:),V(2,:),V(3,:),'*k');
hold off;

[~,I1]=min(Col1(:))
[~,I2]=min(Col2(:))