function plotImageWithRois(Image, Rois, sesname, roi_to_highlight)

if nargin < 4
    roi_to_highlight = [];
end
numrois = size(Rois, 2);

% figure
% imagesc(Image)
Image = mat2gray(Image);
Image = adapthisteq(Image);

imshow(Image)
hold on
for iRoi = 3:numrois
    [~, xx, yy] = roipoly(Image, Rois{iRoi}.mnCoordinates(:, 1), Rois{iRoi}.mnCoordinates(:, 2));
    if ~isempty(roi_to_highlight) && iRoi == roi_to_highlight
        plot(xx, yy, 'LineWidth', 2, 'LineStyle','-.', 'Color', 'r')
    else
        plot(xx, yy, 'LineWidth', 0.8, 'Color', 'g')
    end
    text(mean(xx), mean(yy), num2str(iRoi), 'Color', 'r', 'HorizontalAlignment','center')
end
yticks([])
xticks([])
set(gcf, 'Resize', 'on');
% set(gcf, 'WindowState', 'maximized')
hold off
title(sesname, 'Interpreter','none')