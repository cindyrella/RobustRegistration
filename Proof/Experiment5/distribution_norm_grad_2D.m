%% Compute the norm of gradient of L1 objective function
cmap = jet(numel(0.5:0.5:(2*pi)));
%% Dirty Case
figure()
hold on;
for theta = 0.5:0.5:(2*pi)
    n = 1000;
    gradDirty = zeros(1,n);
    for i = 1:n
        m = 4;
        X = randn(m,2);
        Y = randn(m,2);
        
        gradDirty(i) = norm_grad(X,Y,theta);
    end
    h1 = histogram(gradDirty,'DisplayName',num2str(theta),'FaceColor',cmap(floor(theta*2),:));
end
legend('show')
title(strcat('Distribution gradient for Noisy part small sample'));
saveas(gcf,'Grad_norm_Noisy_small.png');
hold off;
%% Clean Case
figure()
hold on;
for theta = 0.5:0.5:(2*pi)
    n = 1000;
    gradClean = zeros(1,n);
    for i = 1:n
        m = 4;
        X_clean = randn(m,2);
        gradClean(i) = norm_grad(X_clean,X_clean,theta);
    end
    h2 = histogram(gradClean,'DisplayName',num2str(theta),'FaceColor',cmap(floor(theta*2),:));
end
legend('show')
title(strcat('Distribution gradient for Clean part small sample'));
saveas(gcf,'Grad_norm_Clean_small.png');
hold off;