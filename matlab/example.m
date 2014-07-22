img = imread('test.png'); % Read example input image

% ****DON'T FORGET TO CAST TYPE TO DOUBLE****
% This is needed for the ContourTree executable
img = double(img);

[ct_pruning_off ct_img_pruning_off] = ContourTree(img);        % Computes contour tree with no pruning
[ct_pruning_on  ct_img_pruning_on]  = ContourTree(img,0.1,32); % Prune with intensity threshold 0.1 and area/volume threshold 30

% Plot results;
subplot(2,3,1)
imagesc(img);
title('Input image');
axis square;
axis off;

subplot(2,3,2)
imagesc(img);
[x,y] = ind2sub(size(img),ct_pruning_off(:,1));
hold on;
plot(y,x,'or','MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 2);
for i=1:size(ct_pruning_off,1)
    [x1,y1] = ind2sub(size(img),ct_pruning_off(i,1));
    [x2,y2] = ind2sub(size(img),ct_pruning_off(i,2));
    plot([y1 y2], [x1 x2], 'k', 'LineWidth', 2);
end
axis square;
axis off;

subplot(2,3,3)
imagesc(img);
[x,y] = ind2sub(size(img),ct_pruning_on(:,1));
hold on;
plot(y,x,'or','MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 2);
for i=1:size(ct_pruning_on,1)
    [x1,y1] = ind2sub(size(img),ct_pruning_on(i,1));
    [x2,y2] = ind2sub(size(img),ct_pruning_on(i,2));
    plot([y1 y2], [x1 x2], 'k', 'LineWidth', 2);
end
axis square;
axis off;


subplot(2,3,5)
imagesc(ct_img_pruning_off);
title('Contour tree with no pruning');
axis square;
axis off;

subplot(2,3,6)
imagesc(ct_img_pruning_on);
title('Contour tree with pruning. Intensity thresh:0.1,Area/volume thresh:32');
axis square;
axis off;

% Plot contour trees
ContourTreePlot(img,ct_pruning_off,ct_img_pruning_off);
ContourTreePlot(img,ct_pruning_on,ct_img_pruning_on);

