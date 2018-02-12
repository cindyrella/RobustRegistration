clear all
rng(2,'twister');
%% Define problem parameters
m = 10;
psize = 30;
dreal = 3;
x = randn(dreal,psize);
%eps = 1E-8;
nstep = 20;
stepsize = 0.05;
T = zeros(dreal,m);
iter_opt = 200;
cmap = jet(m);

%% Iterate over different noise level
err = zeros(nstep,4);
for l=1:nstep
    %% Simulate patches 
    for i=1:m
        t              = 1*randn(dreal,1)*ones(1,psize); 
        patch(i).idx   = 1:psize;
        noise          = ones(dreal,1)*binornd(1,stepsize*(l-1),1,psize);
        noise          = 1*randn(dreal,psize).*noise;
        patch(i).coord = x+noise+t;
        for j=(i+1):m
            weight(i,j).w = ones(1,psize);
        end
        T(:,i)= t(:,1);
    end
    T = T(:,2:end)-T(:,1);
    
%% Find the LS solution
    [G,Linv,B] = LSSDP(patch,weight,dreal);
    
%% Iterate over different liftings
    ET = zeros(iter_opt+1,m);
    figure()
    for kappa = 1:m
        dlift = dreal*kappa;
        [~,R0] = find_solution_LS(Linv,B,dlift,G);
        
        %% Solve the Riemannian problem
        [eR,eT] = solve_opt_man(patch,dreal,dlift,R0,T,iter_opt);
        semilogy(eR,'DisplayName',strcat('dl =',num2str(dlift)),'Color',cmap(kappa,:));
        hold on;
        ET(:,kappa) = eT';
    end
    xlabel('Iter');
    ylabel('||RR^t-R_{op}R_{op}^t||_F');
    legend('show');
    title(strcat('Error in R with corruption of ',num2str(100*stepsize*(l-1)),'%'));
    saveas(gcf,strcat('eR',num2str(l),'.png'));
    hold off
    figure()
    for kappa = 1:m
        dlift = dreal*kappa;
        semilogy(ET(:,kappa),'DisplayName',strcat('dl =',num2str(dlift)),'Color',cmap(kappa,:));
        hold on;
    end
    xlabel('Iter');
    ylabel('||T-T_{op}||_F');
    legend('show');
    title(strcat('Error in T with corruption of ',num2str(100*stepsize*(l-1)),'%'));
    saveas(gcf,strcat('eT',num2str(l),'.png'));
    hold off
    close all
end
