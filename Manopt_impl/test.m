clear all
%% Define different sources of error
m = 10;
psize = 30;
dreal = 3;
dlift = 6;
x = randn(dreal,psize);
eps = 1E-8;
nstep = 20;
stepsize = 0.05;
T = zeros(dreal,m);

%iterate over different noise level
Err = cell(1,nstep);
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
    [Tst,R] = find_solution_LS(Linv,B,dlift,G);
    niter = 30;
    err(1,1) = norm(R*R'-repmat(eye(dreal),m,m),'fro');
    err(1,2) = norm([T;zeros(dlift-dreal,m-1)]-Tst,'fro');
    
%% Solve the Riemannian problem

% Create manifold
manifold = stiefelstackedfactory(m, dreal, dlift);
problem.M = manifold;

% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(R,store) cost_function(R,patch,dreal,store);
problem.egrad = @(R,store) euc_grad(R,patch,dreal,store);

% Solve.
[x, xcost, info, options] = trustregions(problem,R);
 
% Display some statistics.
figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');
end
