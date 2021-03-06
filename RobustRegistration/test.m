clear all
rng(2,'twister');
%% Define problem parameters
m = 10;
psize = 100;
dreal = 3;

tol = 1E-14;
nstep = 15;
stepsize = 0.05;

ncorr = 7;
corrsize = 0.1;

niter = 30;
tol2 = 1E-8;
cmap = jet(m);

%% Create patches
x = randn(dreal,psize);
Rop = zeros(dreal,dreal * m);
Top = zeros(dreal,m);
for i=1:m
    t   = 1*randn(dreal,1)*ones(1,psize); 
    while true
        [Q,~] = qr(randn(dreal));
        if (det(Q))>tol
            break;
        end
    end
    if (det(Q)+1)<tol
        Q(:,1) = -1*Q(:,1);
    end
        
    patch_op(i).coord = Q'*(x - t);
    for j=(i+1):m
        weight0(i,j).w = ones(1,psize);
    end
    Top(:,i)= t(:,1);
    Rop(:,(dreal*(i-1)+1):(dreal*i)) = Q;
end
Top = Rop(:,1:dreal)'*(Top(:,2:end) - Top(:,1));
Gop = Rop'*Rop;

%% Iterate over different correspondences
for corr = 1:ncorr 
    for i = 1:m
        patch_op(i).idx   = nonzeros((1:psize).*binornd(1,1-corrsize*(corr-1),1,psize));
    end

    %% Iterate over different noise level
    err = zeros(nstep,4);
    for l=1:nstep
        %% Polute patches 
        for i = 1:m
            patch(i).idx   = patch_op(i).idx;
            x              = patch_op(i).coord(:,patch_op(i).idx);
            xsize          = numel(patch_op(i).idx); 
            noise          = ones(dreal,1)*binornd(1,stepsize*(l-1),1,xsize);
            noise          = 1*randn(dreal,xsize).*noise;
            patch(i).coord = x + noise;
        end

    %% Find the LS solution
        [G,Linv,B] = LSSDP(patch,weight0,dreal);

    %% Iterate multiple times
    
        eR = zeros(niter+1,1);
        eT = zeros(niter+1,1);
        [weight, patch,Tst,Gst] = updateweight(patch,Linv,B,m,dreal,eps,weight0,G);
        eR(1) = norm(Gst-Gop,'fro');
        eT(1) = norm(Tst-Top,'fro');
    
        tol3 = tol2;
        for k = 1:niter
            [G,Linv,B] = LSSDP(patch,weight,dreal);
            [weight,patch,Tst,Gst]     = updateweight(patch,Linv,B,m,dreal,tol3,weight,G);
            tol3 = min(tol3/3,1E-16);
            eR(k+1) = norm(Gst-Gop,'fro');
            eT(k+1) = norm(Tst-Top,'fro');
        end
        figure()
        semilogy(eR,'DisplayName','ILS');
        xlabel('Iter');
        ylabel('||RR^t-R_{op}R_{op}^t||_F');
        axis([ 0 35 1E-9 10]);
        legend('show');
        title(strcat('R Error, corruption = ',num2str(100*stepsize*(l-1)),'%, correspondance = ',num2str(100*(1-corrsize*(corr-1))),'%'));
        saveas(gcf,strcat('eR',num2str(l),'_',num2str(corr),'.png'));
        hold off
        figure()
        semilogy(eT,'DisplayName','ILS');
        xlabel('Iter');
        ylabel('||T-T_{op}||_F');
        legend('show');
        axis([ 0 35 1E-9 10]);
        title(strcat('T Error, corruption = ',num2str(100*stepsize*(l-1)),'%, correspondance = ',num2str(100*(1-corrsize*(corr-1))),'%'));
        saveas(gcf,strcat('eT',num2str(l),'_',num2str(corr),'.png'));
        hold off
        close all
    end
end
