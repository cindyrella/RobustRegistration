clear all
rng(2,'twister');
addpath ../
%% Problem variables
m = 10;

%% Read patches
Rop = zeros(dreal,dreal * m);
Top = zeros(dreal,m);
for i = 1:m
    load('../Bunny/view_10.mat',{'index','Q','t','V'})
    patch(i).coord = V;
    for j = (i+1):m
        weight0(i,j).w = ones(1,psize);
    end
    Top(:,i) = t;
    Rop(:,(dreal*(i-1)+1):(dreal*i)) = Q;
end
Top = Rop(:,1:dreal)'*(Top(:,2:end) - Top(:,1));
Gop = Rop'*Rop;


%% Find the LS solution
[G,Linv,B] = LSSDP(patch,weight0,dreal);

%% Iterate over different liftings
ET = zeros(iter_opt+1,m+1);
figure()
%% No lifting
[~,R0] = find_solution_LS(Linv,B,dreal,G);

% Solve the Riemannian problem
[eR,eT] = solve_opt_man(patch,dreal,dreal,R0,Gop,Top,iter_opt);
semilogy(eR,'DisplayName',strcat('No lifting'),'Color','k');
hold on;
ET(1:numel(eT),1) = eT';

%% Lifting

for kappa = 2:m
    dlift = dreal*kappa;
    [~,R0] = find_solution_LS(Linv,B,dlift,G);

    % Solve the Riemannian problem
    [eR,eT] = solve_opt_man(patch,dreal,dlift,R0,Gop,Top,iter_opt);
    semilogy(eR,'DisplayName',strcat('dl =',num2str(dlift)),'Color',cmap(kappa,:));
    ET(1:numel(eT),kappa) = eT';
end

%% ILS
eR = zeros(niter+1,1);
eT = zeros(niter+1,1);
ILS_iter = zeros(niter+1,1);

[weight, patch,Tst,Gst] = updateweight(patch,Linv,B,m,dreal,eps,weight0,G);
eR(1) = norm(Gst-Gop,'fro');
eT(1) = norm(Tst-Top,'fro');
ILS_iter(1) = 1;

tol3 = tol2;
for k = 1:niter
    [G,Linv,B,myiter] = LSSDP(patch,weight,dreal);
    [weight,patch,Tst,Gst]     = updateweight(patch,Linv,B,m,dreal,tol3,weight,G);
    tol3 = max(tol3/3,1E-16);
    eR(k+1) = norm(Gst-Gop,'fro');
    eT(k+1) = norm(Tst-Top,'fro');
    ILS_iter(k+1) = myiter + ILS_iter(k);
end
semilogy(ILS_iter,eR,'DisplayName','ILS','Color','k','LineStyle','--');
ET(1:numel(eT),m+1) = eT;
%% Plot error in R
xlabel('Iter');
ylabel('||RR^t-R_{op}R_{op}^t||_F');
axis([ 0 inf 1E-8 10]);
legend('show');
title(strcat('R Error, corruption = ',num2str(100*stepsize*(l-1)),'%, correspondance = ',num2str(100*(1-corrsize*(corr-1))),'%'));
saveas(gcf,strcat('eR',num2str(l),'_',num2str(corr),'.png'));
hold off
%% Plot error in T
figure()
semilogy(ET(:,1),'DisplayName',strcat('No lifting'),'Color','k');
 hold on;
for kappa = 2:m
    dlift = dreal*kappa;
    semilogy(ET(:,kappa),'DisplayName',strcat('dl =',num2str(dlift)),'Color',cmap(kappa,:));
    hold on;
end

semilogy(ILS_iter,ET(1:(1+niter),m+1),'DisplayName','ILS','Color','k','LineStyle','--');

xlabel('Iter');
ylabel('||T-T_{op}||_F');
legend('show');
axis([ 0 inf 1E-8 10]);
title(strcat('T Error, corruption = ',num2str(100*stepsize*(l-1)),'%, correspondance = ',num2str(100*(1-corrsize*(corr-1))),'%'));
saveas(gcf,strcat('eT',num2str(l),'_',num2str(corr),'.png'));
hold off
close all

