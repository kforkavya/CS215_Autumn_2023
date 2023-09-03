% UpdateMean
% Updates the mean of a dataset with a new data point.

function newMean = UpdateMean(OldMean, NewDataValue, n)
    % Calculate the updated mean using the formula
    newMean = ((n * OldMean) + NewDataValue) / (n + 1);
end

% UpdateMedian
% Updates the median of a dataset with a new data point.

function newMedian = UpdateMedian(oldMedian, NewDataValue, A, n)
    % Insert the new data point into the sorted dataset
    sortedA = sort([A, NewDataValue]);
    
    % Check if the total number of data points is odd or even
    if mod(n + 1, 2) == 0
        % Calculate the average of the two middle values
        newMedian = (sortedA((n+1) / 2) + sortedA(((n+1) / 2) + 1)) / 2;
    else
        % Select the middle value
        newMedian = sortedA((n + 2) / 2);
    end
end

% UpdateStd
% Updates the standard deviation of a dataset with a new data point.

function newStd = UpdateStd(OldMean, OldStd, NewMean, NewDataValue, n)
    % Calculate the updated standard deviation using the formula
    newStd = sqrt(((n-1)*(OldStd)^2 + n*(OldMean)^2 + NewDataValue^2 - (n+1)*(NewMean)^2) / n);
end