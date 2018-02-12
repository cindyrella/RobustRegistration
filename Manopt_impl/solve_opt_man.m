function [Rst,Tst] = solve_opt_man(patch,dreal,dlift,R0)
m = numel(patch);
% Create manifold
manifold = stiefelstackedfactory(m, dreal, dlift);
problem.M = manifold;

% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(R,store) cost_function(R,patch,dreal,store);
problem.egrad = @(R,store) euc_grad(R,patch,dreal,store);
problem.ehess = @(R,eta,store) euc_hess(R,eta,patch,dreal,store);
% Solve.
options.maxiter = 1000;
options.debug = 0;
options.Delta_bar = 0.1;
[Rlift, ~, ~, ~] = trustregions(problem,R0,options);

%%Get projection of lift
[U,sigma] = eigs(Rlift*Rlift',dreal);
O         = U*sqrt(sigma);
[~,O]     = qr(O');
 O        = diag(sign(diag(O)))* O;
 Rst = O';
Tst = find_best_T(Rst,patch,dreal);
Tst = -1*(Tst(:,2:end)-Tst(:,1));
end 