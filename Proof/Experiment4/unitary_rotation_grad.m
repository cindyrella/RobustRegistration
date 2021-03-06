%% Analize unitary rotation

%% Generate data
p = 0.5;
n = 1000;
X = randn(n,3);
P = binornd(1,p,n,1);
Y = randn(n,3).*P+X.*(1-P);
tol = 1E-15;

%% Find principal directions
[~,S,V] = svd(X,0);

%% Precompute entries
% Difference
XmY = diag((X-Y)*(X-Y)');
YcX = cross(Y,X);
XdY = diag(X*Y');

%% Define directions
V1 = sum(YcX)';
V1 = V1/norm(V1);
V2 = cross(cross(V1,V(:,1)),V1);
V2 = V2/norm(V2);
%% Project in directions
V1YcX = YcX * V1;
V2YcX = YcX * V2;

V1X = X * V1;
V2X = X * V2;

V1Y = Y * V1;
V2Y = Y * V2;

%% Objective function
f = @(a,b) XmY - 2*sin(sqrt(a^2+b^2))/max(sqrt(a^2+b^2),tol)*(a*V1YcX + b*V2YcX) + 2*(1-cos(sqrt(a^2+b^2)))*(XdY - 1./max(a^2+b^2,tol)*(a^2*V1X.*V1Y+a*b*(V1Y.*V2X+V1X.*V2Y)+b^2*V2Y.*V2X));

%% Gradient w
fw     = @(a,b) -2*((1-cos(sqrt(a^2+b^2)))/max(sqrt(a^2+b^2),tol)*(Y.*(a*V1X+b*V2X)+X.*(a*V1Y+b*V2Y))+sin(sqrt(a^2+b^2))*YcX);
ftheta = @(a,b) - 2*cos(sqrt(a^2+b^2))/max(sqrt(a^2+b^2),tol)*(a*V1YcX + b*V2YcX) + 2*sin(sqrt(a^2+b^2))*(XdY - 1./max(a^2+b^2,tol)*(a^2*V1X.*V1Y+a*b*(V1Y.*V2X+V1X.*V2Y)+b^2*V2Y.*V2X));

projw  = @(a,b) eye(3) - 1/max(a^2+b^2,tol)*(a^2*(V1*V1')+b^2*(V2*V2')+a*b*(V1*V2'+V2*V1')); 

%% Generate angles


m = 100;
lA = pi/5;
lB = pi/5;
 
vA = linspace(-lA,lA,2*m);
vB = linspace(-lB,lB,2*m);
[B, A] = meshgrid(vB,vA);

%% Add differences

Col1 = zeros(numel(vA),numel(vB));
Col2 = zeros(numel(vA),numel(vB));

Norm_Grad1 = zeros(numel(vA),numel(vB));
Norm_Grad2 = zeros(numel(vA),numel(vB));

Norm_dpsi1 = zeros(numel(vA),numel(vB));
Norm_dpsi2 = zeros(numel(vA),numel(vB));

for id_a = 1:numel(vA)
    id_a
    for id_b = 1:numel(vB)
        vObj = f(vA(id_a),vB(id_b));
        Col1(id_a,id_b) = sum(sqrt(vObj));
        Col2(id_a,id_b) = sum(vObj);
        
        vecw = fw(vA(id_a),vB(id_b));
        Norm_Grad1(id_a,id_b) = norm(sum(vecw./max(sqrt(vObj),tol))*projw(vA(id_a),vB(id_b)));
        Norm_Grad2(id_a,id_b) = norm(sum(vecw)*projw(vA(id_a),vB(id_b)));

        vectheta = ftheta(vA(id_a),vB(id_b));
        Norm_dpsi1(id_a,id_b) = abs(sum(vectheta./max(sqrt(vObj),tol)));
        Norm_dpsi2(id_a,id_b) = abs(sum(vectheta));
    end
end
%% Plot Norms
figure()
subplot(2,3,1);
imagesc([-lA,lA],[-lB,lB],Col1);
hold on;
plot(pi*cos(linspace(0,2*pi,100)),pi*sin(linspace(0,2*pi,100)));
hold off;

subplot(2,3,4);
imagesc([-lA,lA],[-lB,lB],Col2);
hold on;
plot(pi*cos(linspace(0,2*pi,100)),pi*sin(linspace(0,2*pi,100)));
hold off;

subplot(2,3,2);
imagesc([-lA,lA],[-lB,lB],Norm_dpsi1);
hold on;
plot(pi*cos(linspace(0,2*pi,100)),pi*sin(linspace(0,2*pi,100)));
hold off;

subplot(2,3,5);
imagesc([-lA,lA],[-lB,lB],Norm_dpsi2);
hold on;
plot(pi*cos(linspace(0,2*pi,100)),pi*sin(linspace(0,2*pi,100)));
hold off;

subplot(2,3,3);
imagesc([-lA,lA],[-lB,lB],(Norm_Grad1));
hold on;
plot(pi*cos(linspace(0,2*pi,100)),pi*sin(linspace(0,2*pi,100)));
hold off;

subplot(2,3,6);
imagesc([-lA,lA],[-lB,lB],(Norm_Grad2));
hold on;
plot(pi*cos(linspace(0,2*pi,100)),pi*sin(linspace(0,2*pi,100)));
hold off;