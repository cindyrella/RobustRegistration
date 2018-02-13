M = 10;
psize = 30;
d = 3;
x = randn(d,psize);
eps = 1E-8;
nstep = 20;
stepsize = 0.05;
T = zeros(d,M);

%iterate over different noise level
Err = cell(1,nstep);
for l=1:nstep
    
    %Simulate patches 
    for i=1:M
        t              = 1*randn(d,1)*ones(1,psize); 
        patch(i).idx   = 1:psize;
        noise          = ones(d,1)*binornd(1,stepsize*(l-1),1,psize);
        noise          = 1*randn(d,psize).*noise;
        patch(i).coord = x+noise+t;
        for j=i+1:M
            weight(i,j).w = ones(1,psize);
        end
        T(:,i)= t(:,1);
    end
    T = T(:,2:end)-T(:,1);
    
    
    %Start IRLS
    [G,Linv,B] = LSSDP(patch,weight,d);
    niter = 30;
    err = zeros(niter,2);
    
    [weight, patch,Tst] = updateweight(patch,Linv,B,M,d,eps,weight,G);
    err(1,1) = norm(G-repmat(eye(d),M,M),'fro');
    err(1,2) = norm(T-Tst,'fro');
    
    oldG = zeros(size(G));
    for k=1:niter
        if norm(oldG-G,'fro')<10^-5
            break
        end
        oldG       = G;
        [G,Linv,B] = LSSDP(patch,weight,d);
        [weight,patch,Tst]     = updateweight(patch,Linv,B,M,d,eps,weight,G);
        eps = eps/3;
        err(k+1,1) = norm(G-repmat(eye(d),M,M),'fro');
        err(k+1,2) = norm(T-Tst,'fro');
    end
Err{l}=err;
end
