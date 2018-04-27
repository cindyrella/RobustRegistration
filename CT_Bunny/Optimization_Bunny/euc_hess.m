function [sol,store] = euc_hess(R,eta,patch,d,store)
if (~isfield(store, 'hess'))
    [~,store] = euc_grad(R,patch,d,store);
end
hess = store.hess;
sol = trace(eta'*hess*R);
end