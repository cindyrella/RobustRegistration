function [G,store] = euc_grad(X,patch,d,store)
tol = 1E-16;
m = numel(patch);

R = X.R;
T = X.T;

if (~isfield(store, 'Rgrad'))||(~isfield(store, 'Tgrad'))
% Recover the transformations 
if (~isfield(store, 'transf_x'))
    transf_x = cells(1,m);
    for i = 1:m
        xi   = patch(i).coord;
        Ri = X.R((d*(i-1)+1):(d*i),:)';
        ti = X.T(:,i);
        transf_x{i} = Ri * xi - ti;
    end
    store.transf_x = transf_x;
end
transf_x = store.transf_x;

%Compute Grads
Rgrad = zeros(size(R));
Tgrad = zeros(size(T));

for i = 1:m
    idxi = patch(i).idx;
    xi   = patch(i).coord;
    ti = T(:,i);
    for j = (i+1):m
        idxj = patch(j).idx;
        xj   = patch(j).coord;
        [~,intidxi,intidxj] = intersect(idxi,idxj);
        tj = T(:,j);
        
        fxi = transf_x{i}(intidxi);
        xik = xi(intidxi);
        fxj = transf_x{j}(intidxj);
        xjk = xj(intidxj);
        
        errors = fxi - fxj;
        w = vecnorm(errors  + tol);
        
        %cost = cost + sum(w);
      
        Grad_i = - (fxj +  ti) * diag(1./w) * xik'; 
        Grad_j = - (fxi +  tj) * diag(1./w) * xjk';
        
        Rgrad((d*(i-1)+1):(d*i),:) = Rgrad((d*(i-1)+1):(d*i),:) + Grad_i;
        Rgrad((d*(j-1)+1):(d*j),:) = Rgrad((d*(j-1)+1):(d*j),:) + Grad_j;

        tgrad = errors *(1./w)';
        
        Tgrad(:,i) = Tgrad(:,i) - tgrad;
        Tgrad(:,j) = Tgrad(:,j) + tgrad;
    end 
end
store.Rgrad = Rgrad;
store.Tgrad = Tgrad;
end
G.R = store.Rgrad;
G.T = store.Tgrad;
end