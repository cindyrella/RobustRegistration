function [cost,store] = cost_function(X,patch,d,store)
%Computes the cost function assuming optimal value of t for a fixed R
%R is of the size (d) x (3m), where d>= 3 is the lifting
% and m is the number of sets
%patch is the total set of points. 
%patch(i) corresponds to the ith image, and contains .idx (index of points)
% and .coord cooresponds to the 3d coordinates of the points

if(~isfield(store, 'cost'))

m = numel(patch);

cost = 0;

%First apply the transformation to all the coordinates
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

for i = 1:m
    idxi = patch(i).idx;
    for j = (i+1):m
        idxj = patch(j).idx;
        [~,intidxi,intidxj] = intersect(idxi,idxj);
        errors = transf_x{i}(intidxi) - transf_x{j}(intidxj);
        w = vecnorm(errors);
        cost = cost + sum(w);
    end 
end
store.cost = cost;
end
cost = store.cost;

end