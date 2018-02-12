function stats = mystatsfun(problem, Rlift,patch,dreal,T,stats)
%% Get projection of lift
[U,sigma] = eigs(Rlift*Rlift',dreal);
O         = U*sqrt(sigma);
[~,O]     = qr(O');
 O        = diag(sign(diag(O)))* O;
 Rst = O';
%% Get best T
Tst = find_best_T(Rst,patch,dreal);
Tst = -1*(Tst(:,2:end)-Tst(:,1));
%% Compute error
m = numel(patch);
stats.myeR = norm(Rst*Rst'-repmat(eye(dreal),m,m),'fro');
stats.myeT = norm(T-Tst,'fro');
end