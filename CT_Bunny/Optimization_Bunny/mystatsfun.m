function stats = mystatsfun(problem, Rlift,patch,dreal,Gop,Top, stats)
%% Get projection of lift
[U,sigma] = eigs(Rlift*Rlift',dreal);
O         = U*sqrt(sigma);
[~,O]     = qr(O');
 O        = diag(sign(diag(O)))* O;
 Rst = O';
%% Get best T
Tst = find_best_T(Rst,patch,dreal);
Tst = Tst(:,2:end) - Tst(:,1);
%% Compute error
stats.myeR = norm(Rst*Rst'-Gop,'fro');
stats.myeT = norm(Top-Tst,'fro');
end