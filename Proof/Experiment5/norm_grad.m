function [resul] = norm_grad(X,Y,theta)
A = [0,1;-1,0];
norm_x = diag(X*X');
norm_y = diag(Y*Y');
dot_xy = diag(X*Y');
cross_xy = diag(X*A*Y');

base = sqrt(norm_x + norm_y - 2*sin(theta).*cross_xy - 2*cos(theta).* dot_xy);
num = cross_xy ./ base;
denom = dot_xy ./ base;
resul = -cos(theta)*sum(num)+sin(theta)*sum(denom);
end
