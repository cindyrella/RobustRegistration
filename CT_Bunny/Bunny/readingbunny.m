% nlayers = 361;
% size = 512;
% 
% myCT = zeros(size,size,nlayers);
% 
% for i = 1:nlayers
%     X = fopen(num2str(i));
%     M = fread(X,'uint16');
%     myCT(:,:,i) = reshape(M,size,size);
%     fclose('all');
%     imshow(myCT(:,:,i)>40000);
%     imagesc(myCT(:,:,i));
%     myPlot(i) = getframe(gcf);
%     close all;
% end
% 
% video = VideoWriter('myCt2.avi','Uncompressed AVI');
% open(video);
% writeVideo(video,myPlot);
% close(video);
% [X,Y,Z] = meshgrid(0:0.337891:((size - 1) * 0.337891),0:0.337891:((size - 1) * 0.337891),0:0.5:(0.5 * (nlayers - 1 )));
% 
% index = find(myCT > 40000);
% 
% 
% binCT = [X(index),Y(index),Z(index),myCT(index)];
% 
% figure
% scatter3(binCT(:,1),binCT(:,2),binCT(:,3),1,'filled');


%histogram(myCT);
% 
% mode = find((myCT>12400).*(myCT<12800));
% modeCT = [X(mode),Y(mode),Z(mode),myCT(mode)];
% 
% figure
% scatter3(modeCT(:,1),modeCT(:,2),modeCT(:,3),1,'filled');

%Mode corresponds to the border of the circle it is around 12600 

%  index = find(myCT > 12800);
%  
%  histogram(myCT(index));
%  

%Let us visualize the video again removing the background
% figure
% scatter3(binCT(:,1),binCT(:,2),binCT(:,3),1,binCT(:,4),'filled');

%You want to remove when Y>140 

index = find(binCT(:,1)<140);
index2 = find(binCT(index,4)>62000);
 
binCT2 = binCT(index,:);
binCT2 = binCT2(index2,:);


figure
scatter3(binCT2(:,1),binCT2(:,2),binCT2(:,3),1,binCT2(:,4),'filled');

