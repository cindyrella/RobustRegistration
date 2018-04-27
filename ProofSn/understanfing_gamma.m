%% Understanding gamma
dtheta = 100;

x = linspace(0,2*pi,dtheta);
y = linspace(0,2*pi,dtheta);
z = linspace(0,2*pi,dtheta);
[X,Y,Z]=meshgrid(x,y,z);

G1 = max(max(abs(sin(X./2)),abs(sin(Y./2))),abs(sin(Z./2)));
G2 = max(max(abs(cos(X./2)),abs(cos(Y./2))),abs(cos(Z./2)));
G3 = max(max(abs(sin(X)),abs(sin(Y))),abs(sin(Z)));

G = G1.*G2./G3;