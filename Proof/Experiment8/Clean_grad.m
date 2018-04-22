%% Finding out when grad clean is small

%% Fix probability 
q = 0.01;
N = 1:100;

%% Plot different solutions for this q
figure()
plot(N,ones(numel(N),1),'k','DisplayName','\kappa=1');
hold on
for q = [0.01,0.1]
    for k = [-1,0]
        plot(N,-lambertw(k,-q.^(1./N)*exp(-1)),'DisplayName',strcat('q=',num2str(q),'k=',num2str(k)));
    end
end
xlabel('n(1-p)');
ylabel('\kappa');
legend('show')