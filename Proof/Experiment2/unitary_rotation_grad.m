%% Analize unitary rotation

%% Generate data
p = 0.05;
n = 1000;
X = randn(n,3);
P = binornd(1,p,n,1);
Y = randn(n,3).*P+X.*(1-P);

XY = X'*Y;
% [~,S,V] = svd(X,0);
% U = [V(:,3),V(:,1),V(:,2)];
%% Generate angles
m = 200;
thetavec = linspace(0,2*pi,2*m);
phivec = linspace(0,pi,m);
[ph, th] = meshgrid(phivec,thetavec);

R = ones(size(th)); 

b = R.*sin(ph).*cos(th);
c = R.*sin(ph).*sin(th);
d = R.*cos(ph);

%% Add differences

Grad1 = zeros(numel(thetavec),numel(phivec),3);
Grad2 = zeros(numel(thetavec),numel(phivec),3);

Norm_Grad1 = zeros(numel(thetavec),numel(phivec));
Norm_Grad2 = zeros(numel(thetavec),numel(phivec));

Col1 = zeros(numel(thetavec),numel(phivec));
Col2 = zeros(numel(thetavec),numel(phivec));
max_norm1 = 0;
max_norm2 = 0;
for id_t = 1:numel(thetavec)
    for id_p = 1:numel(phivec)
        vec1 = 0;
         accum1 = 0;
        w = [b(id_t,id_p),c(id_t,id_p),d(id_t,id_p)]';
        for id_X = 1:n
            x = X(id_X,:)';
            y = Y(id_X,:)';
            den = (x-y)'*(x-y) - 2*binornd(1,stepsize*(l-1),1,xsize)
            accum1 = accum1 + sqrt(x'*x - (x'*w)^2);
            vec1 = vec1 + x * ((x'*w)/sqrt(x'*x - (x'*w)^2));
        end
        Col1(id_t,id_p) = accum1;
        Col2(id_t,id_p) = -w'*XX*w;
        %% Project into orthogonal space
        %Grad1(id_t,id_p,:) = vec1;
        Grad1(id_t,id_p,:) = vec1 - (w'*vec1)*w;
        Norm_Grad1(id_t,id_p) = norm(vec1 - (w'*vec1)*w);
        vec2 = XX*w;
        %Grad2(id_t,id_p,:) = vec2;
        Grad2(id_t,id_p,:) = vec2 - (w'*vec2)*w;
        Norm_Grad2(id_t,id_p) = norm(vec2 - (w'*vec2)*w);
        %% Scale
        max_norm1 = max(max_norm1,norm(vec1 - (w'*vec1)*w)); 
        max_norm2 = max(max_norm2,norm(vec2 - (w'*vec2)*w)); 
    end
end
Grad1 = Grad1/max_norm1;
Grad2 = Grad2/max_norm2;
%% Plot Norms
s = surf(b,c,d,Col1);
hold on;
s.EdgeColor = 'none';
scatter3([1,0,0],[0,1,0],[0,0,1],'*k');
quiver3(b,c,d,Grad1(:,:,1),Grad1(:,:,2),Grad1(:,:,3))
hold off;
figure()
s = surf(b,c,d,Col2);
hold on;
s.EdgeColor = 'none';
scatter3([1,0,0],[0,1,0],[0,0,1],'*k');
quiver3(b,c,d,Grad2(:,:,1),Grad2(:,:,2),Grad2(:,:,3))
hold off;

%% Plot Gradient 

figure()
s = surf(b,c,d,(Norm_Grad1));
hold on;
scatter3([1,0,0],[0,1,0],[0,0,1],'*k');
s.EdgeColor = 'none';
hold off;
figure()
s = surf(b,c,d,(Norm_Grad2));
hold on;
scatter3([1,0,0],[0,1,0],[0,0,1],'*k');
s.EdgeColor = 'none';