%% Compute the probability of change of sign of gradient in L2 case
p = linspace(0,1,1000);
theta = linspace(0,pi,1000);

%% Remove extrem cases (avoid dividing by zero)
% p = p(2:(end-1));
% theta = theta(2:(end-1));

[P,Theta] = meshgrid(p,theta);

%% Compute probability different n
epsilon = 0.1;
for k = 0:4
n = 10^k;
miu = 2*n.*(1-P).*sin(Theta);
sigma = sqrt(4*n.*(1-P).*(sin(Theta)).^2+ 2*n.*P);
figure()
Prob = cdf('Normal',(epsilon*n-miu)./sigma,0,1)-cdf('Normal',(-epsilon*n-miu)./sigma,0,1);
s = surf(P,Theta,Prob);
title(strcat('P(|G|<n*eps) n= ',num2str(n),' eps =',num2str(epsilon)));
s.EdgeColor = 'none';
saveas(gcf,strcat('L2_eps_n',num2str(k),'.png'));
end