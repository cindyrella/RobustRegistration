function [Tst,R] = find_solution_LS(Linv,B,d,G)
[U,sigma] = eigs(G,d);
O         = U*sqrt(sigma);
[~,O]     = qr(O');
 O        = diag(sign(diag(O)))* O;
T         = -O*B*Linv;


Tst = T(:,2:end) - T(:,1);
R = O';
end