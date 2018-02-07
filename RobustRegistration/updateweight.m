function [weight,patch] = updateweight(patch,Linv,B,M,d,eps,weight,G)

[U,sigma] = eigs(G,d);
O         = U*sqrt(sigma);
O         = O';
T         = -O*B*Linv;


%Update weight
norm_const = 0;
for i=1:M
    idxi = patch(i).idx;
    xi   = patch(i).coord;
    for j=i+1:M
  
        idxj = patch(j).idx;
        xj   = patch(j).coord;
        [~,intidxi,intidxj] = intersect(idxi,idxj);
        w = weight(i,j).w;
        for k=1:numel(intidxi)
            w(k) = 1/sqrt(norm(O(:,(i-1)*d+1:i*d)*xi(:,intidxi(k))-O(:,(j-1)*d+1:j*d)*xj(:,intidxj(k))+T(:,i)-T(:,j))^2+eps);
        end
        norm_const = norm_const + sum(w);
        weight(i,j).w = w; 
    end 
end

%renormalize
psize = size(patch(1).coord,2);
const = psize*M*(M-1)/2;
for i=1:M
    for j=i+1:M
        weight(i,j).w = const*weight(i,j).w/norm_const;
    end
end

% %shift
% for i=1:M
%    patch(i).coord = patch(i).coord + T(:,i)*ones(1,size(patch(i).coord,2));
% end