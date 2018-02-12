function [cost,store] = cost_function(R,patch,d,store)
%Computes the cost function assuming optimal value of t for a fixed R
%R is of the size (d) x (3m), where d>= 3 is the lifting
% and m is the number of sets
%patch is the total set of points. 
%patch(i) corresponds to the ith image, and contains .idx (index of points)
% and .coord cooresponds to the 3d coordinates of the points

%we will reuse L,Lx to create t, if there are not L,Lx then we will take w
%= 1
tol = 1E-16;
if(~isfield(store, 'Lx'))||(~isfield(store, 'Linv'))||(~isfield(store, 'T'))||(~isfield(store, 'cost'))

m = numel(patch);

T = find_best_T(R,patch,d);


%Compute cost, Lx, L
L = zeros(m);
Lx = zeros(d*m,m);
cost = 0;

for i=1:m
    idxi = patch(i).idx;
    xi   = patch(i).coord;
    Ri = R((d*i-2):(d*i),:)';
    ti = T(:,i);
    for j = (i+1):m
        idxj = patch(j).idx;
        xj   = patch(j).coord;
        [~,intidxi,intidxj] = intersect(idxi,idxj);
        Rj = R((d*j-2):(d*j),:)';
        tj = T(:,j);
        w = zeros(1,numel(intidxi));
        for k = 1:numel(intidxi)
            xik = xi(:,intidxi(k));
            xjk = xj(:,intidxj(k));
            w(k) = norm(Ri * xik + ti - Rj * xjk - tj+tol,2);
        end
        cost = cost + sum(w);
        rec_w = 1./w;
        rec_cost = sum(rec_w);
        
        L(i,i) = L(i,i) + rec_cost; 
        L(j,j) = L(j,j) + rec_cost; 
        
        L(i,j) = L(i,j) - rec_cost;
        L(j,i) = L(i,j);
        
        Sxi = sum(xi(:,intidxi).*rec_w,2);
        Sxj = sum(xj(:,intidxj).*rec_w,2);
        
        Lx((d*(i-1)+1):(d*i),i) = Lx((d*(i-1)+1):(d*i),i) + Sxi;
        Lx((d*(i-1)+1):(d*i),j) = Lx((d*(i-1)+1):(d*i),j) - Sxi;
        
        Lx((d*(j-1)+1):(d*j),i) = Lx((d*(j-1)+1):(d*j),i) - Sxj;
        Lx((d*(j-1)+1):(d*j),j) = Lx((d*(j-1)+1):(d*j),j) + Sxj;
        
    end 
end
Linv = pinv(L);
store.Linv = Linv;
store.Lx = Lx;
store.cost = cost;
store.T = T;
end
cost = store.cost;

end