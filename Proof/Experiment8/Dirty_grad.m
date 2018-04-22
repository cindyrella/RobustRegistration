%% Finding out when grad clean is small

%% Fix probability 
q = 0.01;
N = 1:100;

%% Plot different solutions for this q
figure()
hold on
for q = [0.001,0.01,0.1,0.5]
    for k = -1
        kappa = sqrt( (1 + lambertw(k,-q.^(1./N)*exp(-2)/2)).^2-1); 
        plot(N,kappa,'DisplayName',strcat('q=',num2str(q),'k=',num2str(k)));
    end
end
xlabel('np');
ylabel('\kappa');
legend('show')