function [eR,eT] = solve_opt_man(patch,dreal,dlift,R0,T0,Gop,Top,iter)
m = numel(patch);
%% Create Manifold
elements = struct();
elements.R = stiefelstackedfactory(m, dreal, dlift);
elements.T = euclideanfactory(dlift, m);
manifold = productmanifold(elements);
problem.M = manifold;

%% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(X,store) cost_function(X,patch,dreal,store);
problem.egrad = @(X,store) euc_grad(X,patch,dreal,store);
problem.ehess = @(X,eta,store) euc_hess(X,eta,patch,dreal,store);
%% Solve.
options.maxiter = iter;
options.debug = 0;
options.Delta_bar = 0.1;
options.statsfun = @(problem, X, stats) mystatsfun(problem, X,dreal,Gop,Top, stats);

X0.R = R0;
X0.T = T0;
[~, ~, info, ~] = trustregions(problem,X0,options);

%% Get errors
eR = [info.myeR];
eT = [info.myeT];
end 