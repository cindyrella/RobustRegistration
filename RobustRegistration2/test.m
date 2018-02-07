M = 2;
psize = 30;
d = 3;
x = randn(d,psize);
eps = 0.1;
nstep = 2;
err   = zeros(nstep,2);
stepsize = 1/nstep;

%iterate over different noise level
for l=1:nstep
    
    %Simulate patches 
    for i=1:M
        t              = 1*randn(d,1)*ones(1,psize); 
        %[Q,~]         = qr(randn(d));
        %Q              = randn(d);
        Q              = eye(d);
        patch(i).idx   = 1:psize;
        noise          = ones(d,1)*binornd(1,stepsize*(l),1,psize);
        noise          = 1*randn(d,psize).*noise;
        patch(i).coord = Q*x+noise+t;
        for j=i+1:M
            weight(i,j).w = ones(1,psize);
        end
    end
    
    
    
    %Start IRLS
    [G,Linv,B] = LSSDP(patch,weight,d);
    
    [weight, patch] = updateweight(patch,Linv,B,M,d,eps,weight,G);
    
    niter = 30;
    
    for k=1:niter
        oldG       = G;
        [G,Linv,B] = LSSDP(patch,weight,d);
        if norm(oldG-G,'fro')<10^-5
            break
        end
        [weight,patch]     = updateweight(patch,Linv,B,M,d,eps,weight,G);
        eps = eps/3;
        
    end
    err(l,1) = (l-1)*stepsize;
    err(l,2) = norm(G-repmat(eye(d),M,M),'fro')
    
    
end