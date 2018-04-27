function [grad,store] = euc_grad(R,patch,d,store)
tol = 1E-16;
m = numel(patch);
if (~isfield(store, 'grad'))
%Recover Lx and Linv 
if (~isfield(store, 'Lx')) || (~isfield(store, 'Linv'))|| (~isfield(store, 'T'))
    [~,store] = cost_function(R,patch,d,store);
end
Lx = store.Lx;
Linv = store.Linv;
T = store.T;

%Compute cross product
Xcross = zeros(d*m);
for i=1:m
    idxi = patch(i).idx;
    xi   = patch(i).coord;
    Ri = R((d*(i-1)+1):(d*i),:)';
    ti = T(:,i);
    for j = (i+1):m
        idxj = patch(j).idx;
        xj   = patch(j).coord;
        [~,intidxi,intidxj] = intersect(idxi,idxj);
        Rj = R((d*(j-1)+1):(d*j),:)';
        tj = T(:,j);
        M_ij = zeros(d,d);
        M_ii = zeros(d,d);
        M_jj = zeros(d,d);
        for k = 1:numel(intidxi)
            xik = xi(:,intidxi(k));
            xjk = xj(:,intidxj(k));
            w = norm(Ri * xik + ti - Rj * xjk - tj+tol,2);
            M_ij = M_ij + 1./w * (xik*xjk');
            M_ii = M_ii + 1./w * (xik*xik');
            M_jj = M_jj + 1./w * (xjk*xjk'); 
        end
        Xcross( (d*(i-1)+1):(d*i) , (d*(j-1)+1):(d*j) ) = - M_ij; 
        Xcross( (d*(j-1)+1):(d*j) , (d*(i-1)+1):(d*i) ) = - M_ij'; 
%         Xcross( (d*(i-1)+1):(d*i) , (d*(i-1)+1):(d*i) ) = Xcross( (d*(i-1)+1):(d*i) , (d*(i-1)+1):(d*i) )+ M_ii; 
%         Xcross( (d*(j-1)+1):(d*j) , (d*(j-1)+1):(d*j) ) = Xcross( (d*(j-1)+1):(d*j) , (d*(j-1)+1):(d*j) )+ M_jj; 
        
    end 
end
hess = Xcross - Lx*Linv*Lx';
grad = hess*R;
store.grad = grad;
store.hess = hess;
end
grad = store.grad;

end