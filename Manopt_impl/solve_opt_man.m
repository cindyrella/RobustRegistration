function [eR,eT] = solve_opt_man(patch,dreal,dlift,R0,T,iter)
m = numel(patch);
%% Create Manifold
manifold = stiefelstackedfactory(m, dreal, dlift);
problem.M = manifold;

%% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(R,store) cost_function(R,patch,dreal,store);
problem.egrad = @(R,store) euc_grad(R,patch,dreal,store);
problem.ehess = @(R,eta,store) euc_hess(R,eta,patch,dreal,store);
%% Solve.
options.maxiter = iter;
options.debug = 0;
options.Delta_bar = 0.1;
options.statsfun = @(problem, Rlift, stats) mystatsfun(problem, Rlift,patch,dreal,T, stats);

[~, ~, info, ~] = trustregions(problem,R0,options);

%% Get errors
eR = [info.myeR];
eT = [info.myeT];
end 