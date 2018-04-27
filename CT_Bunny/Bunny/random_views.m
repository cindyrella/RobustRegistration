%% Create random views of bunny

%% Parameters
clear all
% Number of slices
m = 10;

%Percentage of points in each slice
s = 0.2;

%Noise per slice
p = 0.2;

%Variance of noise
sigma_noise = 20;
%Variance of translation
sigma_translation = 20;

%% Load data
load('iso_bunny.mat')

%Total number of points
N = size(binCT2,1);

%Center the coordinates
binCT2 = binCT2 - mean(binCT2);

%% Create each view

%Number of points per slice
n = floor(N * s);

Largeviews = zeros(n * m,4);

for i = 1:m
    %Extract indices
    index = sort(randsample(N,n));
    
    %Create noisy observations
    index_noise = binornd(1,p,n,1);
    
    %Create Rotation and translation
    [Q,~] = qr(randn(3,3));
    
    if det(Q)<-1
        Q = -1*Q;
    end
    
    t = randn(1,3)*sigma_translation;
    
    %Create view
    V = (binCT2(index,1:3)*Q' + t).*(1 - index_noise) + (sigma_noise * randn(n,3)) .* index_noise;
    
    %Save visual image
    figure()
    scatter3(V(:,1),V(:,2),V(:,3),1,index_noise,'filled');
    saveas(gcf,strcat('view_',num2str(i),'.fig'));
    
    %Save all the varibales related
    save(strcat('view_',num2str(i),'.mat'),'V','index','index_noise','Q','t');

    %Save for future plot of data
    Largeviews(((i-1)*n + 1): (i * n),:) = [V,index_noise * m + i];
    
end

%Plot all the data
figure()
scatter3(Largeviews(:,1),Largeviews(:,2),Largeviews(:,3),1,Largeviews(:,4),'filled');
saveas(gcf,'Largeview.fig');


