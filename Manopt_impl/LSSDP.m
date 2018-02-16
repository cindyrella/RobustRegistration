function [G,Linv,B,myiter] = LSSDP(patch,weight,d)

M = numel(patch);
%Construct B
I_M = eye(M);
L=0;
B=0;
C=0;
for i=1:M
    idxi = patch(i).idx;
    xi   = patch(i).coord;
    for j=i+1:M
        idxj = patch(j).idx;
        xj   = patch(j).coord;
        [~,intidxi,intidxj] = intersect(idxi,idxj);
        w = weight(i,j).w;
        for k=1:numel(intidxi)
            D = (kron(I_M(:,i),eye(d))*xi(:,intidxi(k))-kron(I_M(:,j),eye(d))*xj(:,intidxj(k)));
            L = L+w(k)*(I_M(:,i)-I_M(:,j))*(I_M(:,i)-I_M(:,j))';
            B = B+w(k)*(kron(I_M(:,i),eye(d))*xi(:,intidxi(k))-kron(I_M(:,j),eye(d))*xj(:,intidxj(k)))*(I_M(:,i)-I_M(:,j))';
            C = C+w(k)* (D*D');
        end
        
    end 
end

Linv = pinv(L);
C = C-B*Linv*B';
%eig(B*Linv*B')
cvx_begin sdp
    variable G(M*d,M*d) symmetric
    minimize trace(C*G)
    G>=0;
    for i=1:M
        G((i-1)*d+1:i*d,(i-1)*d+1:i*d) == eye(d);
    end
cvx_end
myiter = cvx_slvitr;
end
