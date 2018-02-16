clear all
rng(2,'twister');
%% Define problem parameters
m = 10;
psize = 100;
dreal = 3;

tol = 1E-14;
tol2 = 1E-8;
nstep = 14;
stepsize = 0.05;

ncorr = 7;
corrsize = 0.1;

niter = 20;
iter_opt = 300;
cmap = jet(m);

%% Create patches
x = randn(dreal,psize);
Rop = zeros(dreal,dreal * m);
Top = zeros(dreal,m);
orientation = [ones(1,floor(m/2)),-1*ones(1,m-floor(m/2))];
orientidx = randperm(m);
orientation = orientation(orientidx);
for i=1:m
    t   = 1*randn(dreal,1)*ones(1,psize); 
    %Get random sign of rotation
    while true
        [Q,~] = qr(randn(dreal));
        if abs(det(Q))>tol
            break;
        end
    end
    if abs(det(Q)+orientation(i))<tol
        Q(:,1) = -1*Q(:,1);
    end
        
    patch_op(i).coord = Q'*(x - t);
    for j=(i+1):m
        weight0(i,j).w = ones(1,psize);
    end
    Top(:,i)= t(:,1);
    Rop(:,(dreal*(i-1)+1):(dreal*i)) = Q;
end
Top = Rop(:,1:dreal)'*(Top(:,2:end) - Top(:,1));
Gop = Rop'*Rop;

%% Iterate over different correspondences
for corr = 2:ncorr 
    for i = 1:m
        patch_op(i).idx   = nonzeros((1:psize).*binornd(1,1-corrsize*(corr-1),1,psize));
    end

    %% Iterate over different noise level
    err = zeros(nstep,4);
    for l=1:nstep
        %% Polute patches 
        for i = 1:m
            patch(i).idx   = patch_op(i).idx;
            x              = patch_op(i).coord(:,patch_op(i).idx);
            xsize          = numel(patch_op(i).idx); 
            noise          = ones(dreal,1)*binornd(1,stepsize*(l-1),1,xsize);
            noise          = 1*randn(dreal,xsize).*noise;
            patch(i).coord = x + noise;
        end

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
    end
end
