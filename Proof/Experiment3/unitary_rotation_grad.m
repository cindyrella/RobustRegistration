%% Analize unitary rotation

%% Generate data
p = 0.05;
psi = pi/2;
n = 1000;
X = randn(n,3);
P = binornd(1,p,n,1);
Y = randn(n,3).*P+X.*(1-P);

XY = X'*Y;
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

%% Add differences

Grad1 = zeros(numel(thetavec),numel(phivec),3);
Grad2 = zeros(numel(thetavec),numel(phivec),3);

Norm_Grad1 = zeros(numel(thetavec),numel(phivec));
Norm_Grad2 = zeros(numel(thetavec),numel(phivec));

Norm_dpsi1 = zeros(numel(thetavec),numel(phivec));
Norm_dpsi2 = zeros(numel(thetavec),numel(phivec));

Col1 = zeros(numel(thetavec),numel(phivec));
Col2 = zeros(numel(thetavec),numel(phivec));
max_norm1 = 0;
max_norm2 = 0;
for id_t = 1:numel(thetavec)
    id_t
    for id_p = 1:numel(phivec)
        vec1 = 0;
        accum1 = 0;
        dpsi1 = 0;
        
        vec2 = 0;
        accum2 = 0;
        dpsi2 = 0;
        
        w = U*([b(id_t,id_p),c(id_t,id_p),d(id_t,id_p)]');
        
        for id_X = 1:n
            x = X(id_X,:)';
            y = Y(id_X,:)';
            den = (x-y)'*(x-y) - 2*sin(psi)*w'*cross(y,x) + 2*(1-cos(psi))*(x'*y-(w'*x)*(y'*w));
            accum1 = accum1 + sqrt(den);
            accum2 = accum2 + den;
            
            vec1 = vec1 + 1./sqrt(den)*(-sin(psi)*cross(y,x)-(1-cos(psi))*(x*(y'*w)+y*(x'*w)));
            vec2 = vec2 + 2*(-sin(psi)*cross(y,x)-(1-cos(psi))*(x*(y'*w)+y*(x'*w)));
            
            dpsi1 = dpsi1 + 1./sqrt(den)*(-cos(psi)*w'*cross(y,x) + sin(psi)*(x'*y-(w'*x)*(y'*w)));
            dpsi2 = dpsi2 + 2*(-cos(psi)*w'*cross(y,x) + sin(psi)*(x'*y-(w'*x)*(y'*w)));
        end
        
        Col1(id_t,id_p) = accum1;
        Col2(id_t,id_p) = accum2;
        
        %% Project into orthogonal space
        %Grad1(id_t,id_p,:) = vec1;
        Grad1(id_t,id_p,:) = vec1 - (w'*vec1)*w;
        Norm_Grad1(id_t,id_p) = norm(vec1 - (w'*vec1)*w);
        %Grad2(id_t,id_p,:) = vec2;
        Grad2(id_t,id_p,:) = vec2 - (w'*vec2)*w;
        Norm_Grad2(id_t,id_p) = norm(vec2 - (w'*vec2)*w);
        
        %% Derivative of psi
        Norm_dpsi1(id_t,id_p) = dpsi1;
        Norm_dpsi2(id_t,id_p) = dpsi2;
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
%quiver3(b,c,d,Grad1(:,:,1),Grad1(:,:,2),Grad1(:,:,3))
hold off;
figure()
s = surf(b,c,d,Col2);
hold on;
s.EdgeColor = 'none';
scatter3([1,0,0],[0,1,0],[0,0,1],'*k');
%quiver3(b,c,d,Grad2(:,:,1),Grad2(:,:,2),Grad2(:,:,3))
hold off;

%% Plot Gradient 

figure()
s = surf(b,c,d,log(Norm_Grad1));
hold on;
scatter3([1,0,0],[0,1,0],[0,0,1],'*k');
s.EdgeColor = 'none';
hold off;
figure()
s = surf(b,c,d,log(Norm_Grad2));
hold on;
scatter3([1,0,0],[0,1,0],[0,0,1],'*k');
s.EdgeColor = 'none';

%% Plot Gradient wrt psi

figure()
s = surf(b,c,d,(Norm_dpsi1));
hold on;
scatter3([1,0,0],[0,1,0],[0,0,1],'*k');
s.EdgeColor = 'none';
hold off;
figure()
s = surf(b,c,d,(Norm_dpsi2));
hold on;
scatter3([1,0,0],[0,1,0],[0,0,1],'*k');
s.EdgeColor = 'none';