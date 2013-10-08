function [points] = iterfeature(feature)
    points = [];
    while feature.hasNext()
        feature_element = feature.next();
        try
            if feature_element.hasNext()
                points = horzcat(points, iterfeature(feature_element));
            end
            if length(points) == 36
               disp('') 
            end
        catch me
            %me.message
            points = horzcat(points, feature_element);
        end
    end
    feature.resetIteration();
end
