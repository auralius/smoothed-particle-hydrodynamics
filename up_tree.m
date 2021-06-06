clear all;
close all

% load image and view it
I = imread('up_tree_logo.png');

% extract the coordinate
x = [0 0];
k = 1;
for i = 1:2:size(I,1)
    for j = 1:2:size(I,2)
        if I(i,j) > 0
            x(k,:) = [i j]./200;
            k = k + 1;
        end
    end
end

figure
plot(x(:,1),x(:,2),'b.')
axis equal;

dlmwrite('up_tree_logo.mat', x, 'delimiter', '\t')