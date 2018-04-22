%% Create rotation matrix

rt = pi
Q = [cos(rt), sin(rt);-sin(rt),cos(rt)]';

%% Create data
n = 1000;
X = randn(n,2);
theta = linspace(0,2*pi,1000);
cmap = jet(10);
%% Create rotations
Yop = zeros(size(X));
kappa = 5;
for m = 1:n
    t = mod(m,kappa)/kappa*2*pi;
    R = [cos(t), sin(t);-sin(t),cos(t)];
    Yop(m,:) = X(m,:)*R;
end


%% Test different noise

I_p = 0.05:0.1:1;
figure()
for in_p =1:numel(I_p) 
p = I_p(in_p);
P = binornd(1,p,n,1);
Y = (Yop).*P+X.*(1-P);

f = zeros(1,numel(theta));

for i = 1:numel(theta)
    t = theta(i);
    R = [cos(t), sin(t);-sin(t),cos(t)];
    f(i)=sum(sqrt(diag((X-Y*R)*(X-Y*R)')));
end

plot(theta,f,'Color',cmap(in_p,:),'DisplayName',strcat(num2str(p*100),'%'));
hold on;
end
hold off
legend('show');