function c = calculateROIcentroids(roi_set)

numrois = size(roi_set, 2);

c = nan(numrois, 2);
for iRoi = 3:numrois
    c(iRoi, :) = mean(roi_set{iRoi}.mnCoordinates);
end

end
