function [H,store] = euc_hess(X,eta,patch,d,store)
tol = 1E-16;
m = numel(patch);

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

% Compute Hessian at rho, tau
R = X.R;
T = X.T;
rho = eta.R;
tau = eta.T;

RgradR = zeros(size(R));
TgradR = zeros(size(T));

RgradT = zeros(size(R));
TgradT = zeros(size(T));
for i = 1:m
    idxi = patch(i).idx;
    xi   = patch(i).coord;
    ti = T(:,i);
    rhoi = rho((d*(i-1)+1):(d*i),:);
    taui = tau(:,i);
    
    for j = (i+1):m
        idxj = patch(j).idx;
        xj   = patch(j).coord;
        [~,intidxi,intidxj] = intersect(idxi,idxj);
        tj = T(:,j);
        rhoj = rho((d*(j-1)+1):(d*j),:);
        tauj = tau(:,j);
        
        fxi = transf_x{i}(intidxi);
        xik = xi(intidxi);
        fxj = transf_x{j}(intidxj);
        xjk = xj(intidxj);
        
        errors = fxi - fxj;
        w = vecnorm(errors  + tol);
        
        % First term product RGradR
        wi = sum((fxj +  ti).* (rhoi * xik)); 
        wj = sum((fxi +  tj).* (rhoj * xjk)); 
        wei = taui' * errors; 
        wej = tauj' * errors; 
        
      
        Grad_ii = - (fxj +  ti) * diag(wi./w.^3) * xik'; 
        Grad_ij = - (fxj +  ti) * diag(wj./w.^3) * xik'; 
        Grad_ji = - (fxi +  tj) * diag(wi./w.^3) * xjk'; 
        Grad_jj = - (fxi +  tj) * diag(wj./w.^3) * xjk';
        
        % Second term RgradR
        Grad_ij_2 = - xik * diag(1./w) * xjk' * rho_j; 
        Grad_ji_2 = - xjk * diag(1./w) * xik' * rho_i; 
        
        %Group all the terms
        
        RgradR((d*(i-1)+1):(d*i),:) = RgradR((d*(i-1)+1):(d*i),:) + Grad_ii + Grad_ij + Grad_ij_2;
        RgradR((d*(j-1)+1):(d*j),:) = RgradR((d*(j-1)+1):(d*j),:) + Grad_ji + Grad_jj + Grad_ji_2;


        %% Terms in TgradR
        Grad_ii = errors * (wi./w.^3)'; 
        Grad_ij = errors * (wj./w.^3)'; 
        Grad_ji = - Grad_ii; 
        Grad_jj = - Grad_ij;
        
         % Second term RgradR
        Grad_ii_2 = - rho'* xik * (1./w)'; 
        Grad_ij_2 = - Grad_ii_2; 
        Grad_jj_2 = - rho'* xjk * (1./w)'; 
        Grad_ji_2 = - Grad_jj; 
        
        TgradR(:,i) = TgradR(:,i) + Grad_ii + Grad_ij + Grad_ii_2 + Grad_ij_2;
        TgradR(:,j) = TgradR(:,j) + Grad_ji + Grad_jj + Grad_jj_2 + Grad_ji_2;

         %% Terms in RgradT
        Grad_ii =   (fxj +  ti) * diag(wei./w.^3) * xik'; 
        Grad_ij = - (fxj +  ti) * diag(wej./w.^3) * xik'; 
        Grad_ji =   (fxi +  tj) * diag(wei./w.^3) * xjk'; 
        Grad_jj = - (fxi +  tj) * diag(wej./w.^3) * xjk';
        
         % Second term RgradT
        Grad_ii_2 = - xik * (1./w)'*tau_i'; 
        Grad_ij_2 = - Grad_ii_2; 
        Grad_jj_2 = - xik * (1./w)'*tau_i'; 
        Grad_ji_2 = - Grad_jj; 
        
        RgradT(:,i) = RgradT(:,i) + Grad_ii + Grad_ij + Grad_ii_2 + Grad_ij_2;
        RgradT(:,j) = RgradT(:,j) + Grad_ji + Grad_jj + Grad_jj_2 + Grad_ji_2;
        
        %% Terms in TgradT
        Grad_ii =  errors * (wei./w.^3)' ; 
        Grad_ij = -errors * (wej./w.^3)' ; 
        Grad_ji = -errors * (wei./w.^3)' ; 
        Grad_jj =  errors * (wej./w.^3)' ;
        
         % Second term TgradT
        Grad_ii_2 =   sum(1./w) * tau_i'; 
        Grad_ij_2 = - sum(1./w) * tau_j';
        Grad_jj_2 =   sum(1./w) * tau_j'; 
        Grad_ji_2 = - sum(1./w) * tau_i'; 
        
        TgradT(:,i) = TgradT(:,i) + Grad_ii + Grad_ij + Grad_ii_2 + Grad_ij_2;
        TgradT(:,j) = TgradT(:,j) + Grad_ji + Grad_jj + Grad_jj_2 + Grad_ji_2;
    end 
end    
H.R = RgradR + RgradT;
H.T = TgradR + TgradT;
end