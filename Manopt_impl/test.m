clear all
rng(1,'twister');
%% Define different sources of error
m = 10;
psize = 30;
dreal = 3;
dlift = dreal*1;
x = randn(dreal,psize);
%eps = 1E-8;
nstep = 20;
stepsize = 0.05;
T = zeros(dreal,m);

%iterate over different noise level
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
    [Tst,Rst] = find_solution_LS(Linv,B,dlift,G);
    niter = 30;
    err(l,1) = norm(Rst*Rst'-repmat(eye(dreal),m,m),'fro');
    err(l,2) = norm([T;zeros(dlift-dreal,m-1)]-Tst,'fro');
    
%% Solve the Riemannian problem
[Rst,Tst] = solve_opt_man(patch,dreal,dlift,Rst);

err(l,3) = norm(Rst*Rst'-repmat(eye(dreal),m,m),'fro');
err(l,4) = norm(T-Tst,'fro');
end