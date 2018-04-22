function [theta_op] = approx_angle(X,Y)
tol = 1E-4;
A = [0,1;-1,0];
norm_x = diag(X*X');
norm_y = diag(Y*Y');
dot_xy = diag(X*Y');
cross_xy = diag(X*A*Y');

value_op = sum(sqrt(diag((X-Y)*(X-Y)')));
theta_op = 0;
niter = 100;
for j = -4:4
theta = pi/8*j;
for i = 1:niter
    theta_old = theta;
    [arc_theta,value] = compute_angle(norm_x,norm_y,dot_xy,cross_xy,theta);
    theta = atan(arc_theta);
    if abs(theta - theta_old) < tol
        break;
    end
end
% theta
% value
if  value < value_op
  theta_op = theta;
  value_op = value;
end
theta = pi/8*j;
for i = 1:niter
    theta_old = theta;
    [arc_theta,value] = compute_angle(norm_x,norm_y,dot_xy,cross_xy,theta);
    theta = atan(arc_theta)+pi;
    if abs(theta - theta_old) < tol
        break;
    end
end
% theta
% value
if  value < value_op
  theta_op = theta;
  value_op = value;
end
end
I_theta = linspace(-pi,pi,1000);
for i = 1:1000
   [resul,base] = compute_angle(norm_x,norm_y,dot_xy,cross_xy,I_theta(i));
   Resul(i) = resul;
   Base(i) = base;
end
figure()
plot(I_theta,atan(Resul));
hold on;
plot(I_theta,I_theta);
hold on;
plot(I_theta,I_theta-pi);
hold on;
plot(I_theta,I_theta+pi);
hold off;
figure()
plot(I_theta,Base);
theta_op
end
function [resul,base] = compute_angle(norm_x,norm_y,dot_xy,cross_xy,theta)
base = sqrt(norm_x + norm_y - 2*sin(theta).*cross_xy - 2*cos(theta).* dot_xy);
num = cross_xy ./ base;
denom = dot_xy ./ base;
resul = sum(num)/sum(denom);
base = sum(base);
end
