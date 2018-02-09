%% Kymograph w/ straight longitudinal sections.
% Averages each column of each thresholded image of the root. This array of
% averages becomes a single column of the kymograph. Pass the thresholded 3d
% stack.
%%

function [kymograph, detrendedKymograph] = visualise_wave6(thresholded)

[~, xDim, zDim] = size(thresholded);
kymograph = NaN(xDim, zDim);
detrendedKymograph = NaN(xDim, zDim);

% Average columns of each image and save in kymograph
for i = 1:zDim
    image = thresholded(:, :, i);
    columnMeans = nanmean(image); % average each column of image
    kymograph(:, i) = columnMeans; % save as a column of graph    
end

kymograph(kymograph == 0) = NaN;

% Trim empty rows (or 1 val only) from bottom of kymograph
for k = 1:xDim
    if sum(~isnan(kymograph(end, :))) <= 1
        kymograph(end, :) = [];
        xDim = xDim - 1;
    else
        break
    end
end

% Filter, detrend and normalise by row
for j = 1:xDim
    row = kymograph(j, :);
    
    firstVal = find(~isnan(row), 1, 'first');
    lastVal = find(~isnan(row), 1, 'last');
    row = row(firstVal:lastVal);

    % handle exception where theres zero or one value only.
    if numel(row) < 2
        continue
    end
    
    % interpolate small NaN gaps in the middle due to threshold errors.
    row = interp1gap(row, 3, 'spline');
        
    detrendedVals = amp_detrend(row', 10, -1);
    detrendedVals = detrendedVals';
    [b, a] = butter(3, 0.15); % creates 3rd order butterworth filter
    
    try 
        detrendedVals = filtfilt(b, a, detrendedVals);
    catch % filtfilt doesnt work for <9 data points
        detrendedVals = detrendedVals; 
    end

    % normalise 0-1
    inputMax = max(detrendedVals);
    inputMin = min(detrendedVals);
    detrendedVals = detrendedVals - inputMin;
    detrendedVals = detrendedVals/(inputMax-inputMin);
   
    detrendedKymograph(j, firstVal:lastVal) = detrendedVals;
    kymograph(j, firstVal:lastVal) = row;

end

detrendedKymograph(detrendedKymograph == 0) = NaN;

% trim empty rows again
for k = 1:xDim
    if sum(~isnan(kymograph(end, :))) <= 1
        kymograph(end, :) = [];
        detrendedKymograph(end, :) = [];
        xDim = xDim - 1;
    else
        break
    end
end

end
yTicksA = 1;
yLabelA = -round((hypJ-1) * 0.061, 1);
% 
if hypJ == 1
    yTicksA = [];
    yLabelA = [];
end

% pixel size = 61 um (lumo 4x4)
yTicksB = hypJ:(1/0.061):xDim;
yLabelB = 0:1:rootL*0.061;

%y axis from top of cut root
% yLabel = 0:1:xDim*0.061;
% yTicks = 1:1/0.061:xDim;

% from RT
% yLabel = fliplr(0:1:xDim*0.061);
% yTicks = xDim-(max(yLabel)/0.061):1/0.061:xDim;
% 
% set(gca, 'XTick', xTicks, 'XTickLabel', timeLabel, 'yTick', yTicks,...
% 'YTickLabel', yLabel, 'FontSize', 7)% from root tip
set(gca, 'XTick', xTicks, 'XTickLabel', timeLabel, 'yTick', [yTicksA, yTicksB], 'YTickLabel', [yLabelA, yLabelB], 'FontSize', 7) % from hyp junction
xlabel('Time (h)');
ylabel({'Distance from hypocotyl', 'junction (mm)'});
%colorbar
%ylabel({'Distance from root tip', '(mm)'});
print('kymographs/working', '-dpdf', '-r300', '-painters')

end
