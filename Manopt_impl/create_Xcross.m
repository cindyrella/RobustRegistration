function [Xcross] = create_Xcross(patch)
m = numel(patch);

%Compute cross product
Xcross = zeros(3*m);
for i=1:m
    idxi = patch(i).idx;
    xi   = patch(i).coord;
    for j = (i+1):m
        idxj = patch(j).idx;
        xj   = patch(j).coord;
        [~,intidxi,intidxj] = intersect(idxi,idxj);
        M = zeros(3,3);
        for k = 1:numel(intidxi)
            xik = xi(:,intidxi(k));
            xjk = xj(:,intidxj(k));
            M = M + xik*xjk';
        end
        Xcross( (3*i-2):(3*i) , (3*j-2):(3*j) ) = M; 
        Xcross( (3*j-2):(3*j) , (3*i-2):(3*i) ) = M'; 
        
    end 
end
end