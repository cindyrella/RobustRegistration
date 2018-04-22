%% Finding out when grad clean is small

%% Fix probability 
N = 1:1000;

% %% Plot different solutions for this q
% figure()
% loglog(N,ones(numel(N),1),'k','DisplayName','Permissible');
% hold on
% for q = [0.01,0.1]
%     for p = [0.05,0.1,0.25,0.5]
%         kappa = sqrt( (1 + lambertw(-1,-q.^(1./(p*N))*exp(-2)/2)).^2-1);
%         ratio = -1*p.*kappa./(2*(1-p).*lambertw(0,-exp(-1)*q.^(1./(N*(1-p)))));
%         loglog(N,ratio,'DisplayName',strcat('q=',num2str(q),'p=',num2str(p)));
%     end
% end
% xlabel('n');
% ylabel('sin\theta');
% legend('show')

%% Plot different solutions for this q
figure()
%loglog(N,ones(numel(N),1),'k','DisplayName','Permissible');
hold on
for p = [0.05,0.1,0.25,0.5]
    for q = [0.01,0.1,0.5]
        kappa = sqrt( (1 + lambertw(-1,-q.^(1./(p*N))*exp(-2)/2)).^2-1);
        ratio = -1*p.*kappa./(2*(1-p).*lambertw(0,-exp(-1)*q.^(1./(N*(1-p)))));
        plot(N,180*asin(ratio)/pi,'DisplayName',strcat('p=',num2str(p),'q=',num2str(q)));
    end
end
xlabel('n');
ylabel('\theta (grad)');
legend('show')