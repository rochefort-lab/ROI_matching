function roi_angles = compute_angles(centroids)
    num_centroids = size(centroids, 1);
    roi_angles = zeros(num_centroids, num_centroids);
    
    for i = 1:num_centroids
        for j = 1:num_centroids
            if i ~= j
                roi_angles(i, j) = atan2(centroids(j, 2) - centroids(i, 2), centroids(j, 1) - centroids(i, 1));
            end
        end
    end
end
