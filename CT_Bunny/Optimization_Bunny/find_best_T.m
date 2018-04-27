function [T] = find_best_T(R,patch,d)
% first approx for t
[Lx,Linv] = create_L(patch,d);
%Compute t 
T = -R'*Lx*Linv;

max_iter = 10;
tol = 1E-10;

for iter = 1:max_iter
    Told = T;
    [Lx,Linv] = create_L_wRT(patch,d,R,Told);
    T = -R'*Lx*Linv;
    if norm(T - Told,'fro') < tol
        break;
    end
end

end