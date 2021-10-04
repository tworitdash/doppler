% M-file to make random clusters of points.
% Clean up
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures if you have the Image Processing Toolbox.
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 15;

% Initialize some parameters.
numberOfGaussians = 20; % Number of small Gaussians
bigImageWidth = 500;
bigImageHeight =  500; % square area 0f 500*500

% Make an image of random points.
grayImage = zeros(bigImageHeight, bigImageWidth, 'uint8');
noiseImage = imnoise(grayImage, 'salt & pepper', .1);
subplot(2, 3, 1); 
imshow(noiseImage, []);
axis on;
title('Initially Random Points', 'FontSize', fontSize);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

gaussianWindowWidth = 121; % Width of the window containing the Gaussian.
gaussianStdDev = 25; % Spread of Gaussian within the window.
gaussian1 = fspecial('Gaussian', gaussianWindowWidth, gaussianStdDev);
% Display it in the upper left plot.
subplot(2, 3, 2);
imshow(gaussian1, []);
axis on;
title('Single Gaussian Image', 'FontSize', fontSize);

% Get random coordinates in the big image where
% we will place the upper left corner of the Gaussians.
widthThatWillFit = bigImageWidth - gaussianWindowWidth;
heightThatWillFit = bigImageHeight - gaussianWindowWidth;
% Get all the small Gaussian locations in arrays of their x and y locations (their upper lefts).
smallUpperLeftX = widthThatWillFit * rand(numberOfGaussians, 1);
smallUpperLeftY = heightThatWillFit * rand(numberOfGaussians, 1);
% Initialize an output image to hold many small overlapping Gaussians.
manySmallGaussians = zeros(bigImageHeight, bigImageWidth);
% Place the small Gaussians one by one.
for k = 1 : numberOfGaussians
	% Find the square in the big image where we're going to add a small Gaussian.
	x1 = int16(smallUpperLeftX(k));
	y1 = int16(smallUpperLeftY(k));
	x2 = int16(x1 + gaussianWindowWidth - 1);
	y2 = int16(y1 + gaussianWindowWidth - 1);
	% Add in one small Gaussian to the existing big image.
% 	manySmallGaussians(y1:y2, x1:x2) = manySmallGaussians(y1:y2, x1:x2) + gaussian1;
	% Assign one small Gaussian to the existing big image.
	manySmallGaussians(y1:y2, x1:x2) =  max(manySmallGaussians(y1:y2, x1:x2), gaussian1);
end
% Scale to the 0-1 range.
manySmallGaussians = mat2gray(manySmallGaussians);
% Display it in the lower left plot.
subplot(2, 3, 3);
imshow(manySmallGaussians, []);
axis on;
title('Many Overlapping Gaussians', 'FontSize', fontSize);

% Multiply the big Gaussian mask by the noise.
maskedBySmallGaussians = double(noiseImage) .* manySmallGaussians;
% Display it.
subplot(2, 3, 4);
imshow(maskedBySmallGaussians, []);
axis on;
title('Noise Masked by Clusters', 'FontSize', fontSize);

% Binarize it.
clusterImage = maskedBySmallGaussians > 150;
% Display it.
subplot(2, 3, 5);
imshow(clusterImage, []);
axis on;
title('Binary Clusters', 'FontSize', fontSize);
