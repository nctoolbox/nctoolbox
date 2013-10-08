function [points] = iterstations(stations)
    points = [];
    for i = 1:length(stations)
        station = stations.get(i-1);
        points = horzcat(points, iterfeature(station));
    end
end
