function stats = mystatsfun(problem, X,dreal,Gop,Top, stats)
%% Get values
R = X.R;
T = X.T;

%% Get projection of lift
[U,sigma] = eigs(Rlift*Rlift',dreal);
O         = U*sqrt(sigma);
[~,O]     = qr(O');
 O        = diag(sign(diag(O)))* O;
 Rst = O';
%% Get projection T
Tst = T(1:dreal,:);
Tst = Tst(:,2:end) - Tst(:,1);
%% Compute error
stats.myeR = norm(Rst*Rst'-Gop,'fro');
stats.myeT = norm(Top-Tst,'fro');
end