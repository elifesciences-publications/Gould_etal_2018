%% 
% Analysis code from "Coordination of robust single cell rhythms in  the Arabidopsis % circadian clock via spatial waves of gene expression" (BioRxiv)
% doi https://doi.org/10.1101/208900s 

% Created by Mark Greenwood (U. of Cambridge). 09/02/2018.
%
% 
% This function thresholds the 3d image stack. 
% takes image stack (3d matrix), and the y pixel position of the hypocotyl junction.
% outputs a 3d matrix containing only the thresholded pixels.
%% 

function thresholdedStack=image_filter(imageStack, hypJ)

start = hypJ - 10;
start(start<1) = 1;
if start > 1
    imageStack(:, 1:start, :) = []; 
end

[yDim, xDim, zDim] = size(imageStack);

% output thresholded matrix
thresholdedStack = NaN(yDim, xDim, zDim);

% filter and threshold at image level
for w = 1:zDim
    frame = imageStack(:, :, w);
    frame(frame==0) = NaN;

%   mean threshold
    thresh = nanmean(frame(:));
    frame(frame<thresh) = NaN;
    thresholdedStack(:, :, w) = frame;
    
end

thresholdedStack(thresholdedStack <= 0) = NaN; 

end



