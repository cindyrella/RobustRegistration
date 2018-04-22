%% Testing Chebyshev bounds

p = linspace(0,1,1000);
theta = linspace(0,pi,1000);

%% Remove extreme cases (avoid dividing by zero)
% p = p(2:(end-1));
% theta = theta(2:(end-1));

[P,Theta] = meshgrid(p,theta);

% %% Compute 1/eps_optimal and 1/(1-eps) optimal

% rec_eps = 1+1./(2.*(1./P-1).*(sin(Theta)).^2).^(1/3);
% rec_onemeps = 1+(2.*(1./P-1).*(sin(Theta)).^2).^(1/3);

% %% Compute the total probability (forgetting the factor of 1/n)
% 
% H = P./(sin(Theta)).^2.*rec_onemeps.^2+2.*(1-P).*rec_eps;
% H_1 = P./(sin(Theta)).^2.*rec_onemeps.^2;
% H_2 =  2.*(1-P).*rec_eps;
% H = H .* (sin(Theta)).^2;
% 
% s = surf(P,Theta,H);
% s.EdgeColor = 'none';
% figure()
% s = surf(P,Theta,log(H_1.* (sin(Theta)).^2));
% s.EdgeColor = 'none';
% figure()
% s = surf(P,Theta,log(H_2.* (sin(Theta)).^2));
% s.EdgeColor = 'none';

%% Compute total probability forgetting factor (1/(n*sin(theta)))

H = ((2.*(1-P)).^(1/3).*sin(Theta).^(2/3)+P.^(1/3)).^3;
 s = surf(P,Theta,H);
 s.EdgeColor = 'none';