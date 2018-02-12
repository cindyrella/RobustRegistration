clear all
rng(2,'twister');
%% Define problem parameters
m = 10;
psize = 30;
dreal = 3;

tol = 1E-14;
nstep = 20;
stepsize = 0.05;

iter_opt = 200;
cmap = jet(m);

%% Create patches
x = randn(dreal,psize);
Rop = zeros(dreal,dreal * m);
Top = zeros(dreal,m);
for i=1:m
    t                 = 1*randn(dreal,1)*ones(1,psize); 
    patch_op(i).idx   = 1:psize;
    while true
        [Q,~] = qr(randn(dreal));
        if (det(Q))>tol
            break;
        end
    end
    if (det(Q)+1)<tol
        Q(:,1) = -1*Q(:,1);
    end
        
    patch_op(i).coord = Q'*(x - t);
    for j=(i+1):m
        weight(i,j).w = ones(1,psize);
    end
    Top(:,i)= t(:,1);
    Rop(:,(dreal*(i-1)):(dreal*i)) = Q;
end
Top = Top(:,2:end) - Top(:,1);
Gop = Rop'*Rop;
%% Iterate over different noise level
err = zeros(nstep,4);
for l=1:nstep
    %% Polute patches 
    for i = 1:m
        x              = patch_op(i).coord;
        noise          = ones(dreal,1)*binornd(1,stepsize*(l-1),1,psize);
        noise          = 1*randn(dreal,psize).*noise;
        patch(i).coord = x + noise;
    end
    
%% Find the LS solution
    [G,Linv,B] = LSSDP(patch,weight,dreal);
    
%% Iterate over different liftings
    ET = zeros(iter_opt+1,m);
    figure()
    for kappa = 1:m
        dlift = dreal*kappa;
        [~,R0] = find_solution_LS(Linv,B,dlift,G);
        
        %% Solve the Riemannian problem
        [eR,eT] = solve_opt_man(patch,dreal,dlift,R0,Gop,Top,iter_opt);
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
