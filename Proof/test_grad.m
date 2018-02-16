%% See how the rotation norm works
clear all
rng(2,'twister');
%% Create the points and the noise
d = 3;
m = 100;
x = randn(d,m);
tol = 1E-16;

np = 21;
nstep = 1/(np-1);

%% Create discretization of sphere
%Nice Fibonacci discretization but not easy to plot
% n = 300;
% index = -n:n;
% golden_inv = (sqrt(5)-1)/2;
% lat = arcsin(2.*index./(2.*n+1));
% lon = 2*pi*index*golden_inv;

% Dumb discretization of angles
n = 60;
phi = linspace(0,pi,n);
theta = linspace(0,2*pi,n);
psi = linspace(0,pi,n);

[Theta,Phi, Psi] = meshgrid(theta,phi,psi);


%% First without noise 
% [Q,~] = qr(randn(3));
% if (det(Q))<0
%     Q(:,1)=-1*Q(:,1);
% end
y = x;


Resul = zeros(n,n,n);
for i = 1:n
    for j = 1:n
        for k = 1:n
            % Compute rotation
            %Create axis
            v = [sin(phi(i))*cos(theta(j)),sin(phi(i))*sin(theta(j)),cos(phi(i))];
            A = zeros(3,3);
            A(1,2) = -v(3);
            A(1,3) = v(2);
            A(2,3) = -v(1);
            A = A - triu(A)';
            
            %Compute rotation
            R = expm(A*psi(k));
            %Compute L1 error
            M = 0;
            for l = 1:m
                den = max(norm(x(:,l) - R*y(:,l)),tol);
                M = M + 1./den.*(-x(:,l)*y(:,l)'+(R*y(:,l))*(R'*x(:,l))');
            end
            Resul(i,j,k) = norm(M,'fro');
        end
    end
end
for val = linspace(0,max(Resul(:)),30)
    isosurface(Theta,Phi,Psi,Resul,val);
    hold on;
end
