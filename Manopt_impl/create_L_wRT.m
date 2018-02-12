function [Lx,Linv] = create_L_wRT(patch,d,R,T)
tol = 1E-16;
m = numel(patch);
%Construct L and Lx
L = zeros(m);
Lx = zeros(d*m,m);
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
end