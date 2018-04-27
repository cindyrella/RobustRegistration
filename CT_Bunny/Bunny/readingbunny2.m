 clear all
 nlayers = 361;
 size = 512;
 
 %% Read CT slices
 
 myCT = zeros(size,size,nlayers);
 
for i = 1:nlayers
    X = fopen(num2str(i));
    M = fread(X,'uint16');
    myCT(:,:,i) = reshape(M,size,size);
    fclose('all');
    close all;
end
 
%% Create frame with XYZ coordinates
 [X,Y,Z] = meshgrid(0:0.337891:((size - 1) * 0.337891),0:0.337891:((size - 1) * 0.337891),(0.5 * (nlayers - 1 )):-0.5:0);
 
 %% Take the points with larger magnitude
 index = find(myCT > 62000);
 
binCT = [X(index),Y(index),Z(index),myCT(index)];

%% Remove wall of background

index = find(binCT(:,1)<120);
 
binCT2 = binCT(index,:);

figure
scatter3(binCT2(:,1),binCT2(:,2),binCT2(:,3),1,'filled');
saveas(gcf,'Original_bunny.fig');

save('iso_bunny.mat','binCT2')