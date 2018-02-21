%% See how the rotation norm works
clear all
rng(2,'twister');
%% Create the points and the noise
d = 2;
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
n = 100;
phi = linspace(0,2*pi,n);



%% First without noise 
% [Q,~] = qr(randn(3));
% if (det(Q))<0
%     Q(:,1)=-1*Q(:,1);
% end
y = x;


Resul = zeros(n,4);
for i = 1:n
    R = [cos(phi(i)),sin(phi(i));-sin(phi(i)),cos(phi(i))];
            %Compute L1 error
            sum1 = 0;
            sum2 = 0;
            M1 = 0;
            M2 = 0;
            for l = 1:m
                den = max(norm(x(:,l) - R*y(:,l)),tol);
                M1 = M1 + 1./den.*(-x(:,l)*y(:,l)'+(R*y(:,l))*(R'*x(:,l))');
                sum1 = sum1 + norm(x(:,l) - R*y(:,l));
                sum2 = sum2 + norm(x(:,l) - R*y(:,l))^2;
                M2 = M2 + (-x(:,l)*y(:,l)'+(R*y(:,l))*(R'*x(:,l))');
            end
            Resul(i,1) = sum1;
            Resul(i,2) = sum2;
            Resul(i,3)= norm(M1,'fro');
            Resul(i,4)= norm(M2,'fro');
  
end
plot(phi,Resul(:,1));
hold on;
plot(phi,Resul(:,2));
plot(phi,Resul(:,3));
plot(phi,Resul(:,4));
