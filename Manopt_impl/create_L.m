function [Linv,Lx] = create_L(patch,d)
m = numel(patch);
%Construct L and Lx
L = zeros(m);
Lx = zeros(d*m,m);
for i = 1:m
    idxi = patch(i).idx;
    xi   = patch(i).coord;
    for j = (i+1):m
        idxj = patch(j).idx;
        xj   = patch(j).coord;
        [~,intidxi,intidxj] = intersect(idxi,idxj);
        
        L(i,i) = L(i,i) + numel(intidxi); 
        L(j,j) = L(j,j) + numel(intidxi); 
        
        L(i,j) = L(i,j) - numel(intidxi);
        L(j,i) = L(i,j);
        
        Sxi = sum(xi(:,intidxi),2);
        Sxj = sum(xj(:,intidxj),2);
        
        Lx((d*(i-1)+1):(d*i),i) = Lx((d*(i-1)+1):(d*i),i) + Sxi;
        Lx((d*(i-1)+1):(d*i),j) = Lx((d*(i-1)+1):(d*i),j) - Sxi;
        
        Lx((d*(j-1)+1):(d*j),i) = Lx((d*(j-1)+1):(d*j),i) + Sxj;
        Lx((d*(j-1)+1):(d*j),j) = Lx((d*(j-1)+1):(d*j),j) - Sxj;
        
    end 
end
Linv = pinv(L);
end