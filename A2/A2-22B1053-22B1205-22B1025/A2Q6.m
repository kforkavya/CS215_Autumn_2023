% Load the first image (im) and the second image (im2)
image1 = double(imread('T1.jpg'));
image2 = double(imread('T2.jpg'));

% Define the bin width for histogram calculations
bin_width = 10;

% Define a range of shifts (tx) to explore
shift_values = -10:10;

% Initialize arrays to store correlation coefficients and QMI values
correlation_coefficients = zeros(1, length(shift_values));
qmi_values = zeros(1, length(shift_values));

% Loop through different shift values (tx)
for i = 1:length(shift_values)
    shift_x = shift_values(i);
    
    % Apply the shift to the second image (im2)
    if shift_x < 0
        shifted_image2 = [image2(:, abs(shift_x) + 1:end), zeros(size(image2, 1), abs(shift_x))];
    elseif shift_x > 0
        shifted_image2 = [zeros(size(image2, 1), shift_x), image2(:, 1:end - shift_x)];
    else
        shifted_image2 = image2; % No shift
    end
    
    % Initialize histograms and joint histogram
    num_bins = round(256 / bin_width);
    joint_histogram = zeros(num_bins, num_bins);
    histogram_image1 = zeros(num_bins, 1);
    histogram_shifted_image2 = zeros(num_bins, 1);

    % Calculate joint histogram
    for a = 1:size(image1, 1)
        for j = 1:size(image1, 2)
            bin_image1 = floor(double(image1(a, j)) / bin_width) + 1;
            bin_shifted_image2 = floor(double(shifted_image2(a, j)) / bin_width) + 1;
            joint_histogram(bin_image1, bin_shifted_image2) = joint_histogram(bin_image1, bin_shifted_image2) + 1;
        end
    end

    % Calculating Marginal histograms from Joint Histograms
    for x = 1:size(joint_histogram, 1)
        for y = 1:size(joint_histogram, 2)
            histogram_image1(x) = histogram_image1(x) + joint_histogram(x, y);
            histogram_shifted_image2(y) = histogram_shifted_image2(y) + joint_histogram(x, y);
        end
    end
    
    % Normalize histograms and joint histogram
    joint_histogram = joint_histogram / sum(joint_histogram(:));
    histogram_image1 = histogram_image1 / sum(histogram_image1);
    histogram_shifted_image2 = histogram_shifted_image2 / sum(histogram_shifted_image2);

    % Calculate the correlation coefficient (ρ)
    % Compute the means of each image manually
    mean_image1 = sum(image1(:)) / numel(image1);
    mean_shifted_image2 = sum(shifted_image2(:)) / numel(shifted_image2(:));

    numerator=0;
    denominator_image1=0;
    denominator_shifted_image2=0;
    % Calculate the numerator and denominators manually
    for a = 1:size(image1, 1)
        for j = 1:size(image1, 2)
            diff_image1 = image1(a, j) - mean_image1;
            diff_shifted_image2 = shifted_image2(a, j) - mean_shifted_image2;
            
            numerator = numerator + (diff_image1 * diff_shifted_image2);
            denominator_image1 = denominator_image1 + (diff_image1^2);
            denominator_shifted_image2 = denominator_shifted_image2 + (diff_shifted_image2^2);
        end
    end

    % Finally
    correlation_coefficients(i) = numerator / (sqrt(denominator_image1) * sqrt(denominator_shifted_image2));

    % Calculate the Quadratic Mutual Information (QMI)
    qmi = 0;
    for a = 1:size(joint_histogram, 1)
        for j = 1:size(joint_histogram, 2)
            qmi = qmi + (joint_histogram(a, j) - histogram_image1(a) * histogram_shifted_image2(j))^2;
        end
    end
    qmi_values(i) = qmi;
end

% Create a figure with two subplots to display results
figure;

% Plot 1: Correlation Coefficient vs. Shift (Original Image 2)
plot(shift_values, correlation_coefficients, '-o');
xlabel('Shift (tx pixels)');
ylabel('Correlation Coefficient');
title('Correlation Coefficient vs. Shift');

figure;
% Plot 2: QMI vs. Shift (Original Image 2)
plot(shift_values, qmi_values, '-o');
xlabel('Shift (tx pixels)');
ylabel('Quadratic Mutual Information (QMI)');
title('QMI vs. Shift');

% Same process for Negative of Image 1
image2 = 255 - image1;

% Initialize arrays to store correlation coefficients and QMI values
correlation_coefficients = zeros(1, length(shift_values));
qmi_values = zeros(1, length(shift_values));

% Loop through different shift values (tx)
for i = 1:length(shift_values)
    shift_x = shift_values(i);
    
    % Apply the shift to the second image (im2)
    if shift_x < 0
        shifted_image2 = [image2(:, abs(shift_x) + 1:end), zeros(size(image2, 1), abs(shift_x))];
    elseif shift_x > 0
        shifted_image2 = [zeros(size(image2, 1), shift_x), image2(:, 1:end - shift_x)];
    else
        shifted_image2 = image2; % No shift
    end
    
    % Initialize histograms and joint histogram
    num_bins = round(256 / bin_width);
    joint_histogram = zeros(num_bins, num_bins);
    histogram_image1 = zeros(num_bins, 1);
    histogram_shifted_image2 = zeros(num_bins, 1);

    % Calculate joint histogram
    for a = 1:size(image1, 1)
        for j = 1:size(image1, 2)
            bin_image1 = floor(double(image1(a, j)) / bin_width) + 1;
            bin_shifted_image2 = floor(double(shifted_image2(a, j)) / bin_width) + 1;
            joint_histogram(bin_image1, bin_shifted_image2) = joint_histogram(bin_image1, bin_shifted_image2) + 1;
        end
    end

    % Calculating Marginal histograms from Joint Histograms
    for x = 1:size(joint_histogram, 1)
        for y = 1:size(joint_histogram, 2)
            histogram_image1(x) = histogram_image1(x) + joint_histogram(x, y);
            histogram_shifted_image2(y) = histogram_shifted_image2(y) + joint_histogram(x, y);
        end
    end
    
    % Normalize histograms and joint histogram
    joint_histogram = joint_histogram / sum(joint_histogram(:));
    histogram_image1 = histogram_image1 / sum(histogram_image1);
    histogram_shifted_image2 = histogram_shifted_image2 / sum(histogram_shifted_image2);
    
    numerator=0;
    denominator_image1=0;
    denominator_shifted_image2=0;
    % Calculate the correlation coefficient (ρ)
    % Compute the means of each image manually
    mean_image1 = sum(image1(:)) / numel(image1);
    mean_shifted_image2 = sum(shifted_image2(:)) / numel(shifted_image2(:));

    % Calculate the numerator and denominators manually
    for a = 1:size(image1, 1)
        for j = 1:size(image1, 2)
            diff_image1 = image1(a, j) - mean_image1;
            diff_shifted_image2 = shifted_image2(a, j) - mean_shifted_image2;
            
            numerator = numerator + (diff_image1 * diff_shifted_image2);
            denominator_image1 = denominator_image1 + (diff_image1^2);
            denominator_shifted_image2 = denominator_shifted_image2 + (diff_shifted_image2^2);
        end
    end

    % Finally
    correlation_coefficients(i) = numerator / (sqrt(denominator_image1) * sqrt(denominator_shifted_image2));

    % Calculate the Quadratic Mutual Information (QMI)
    qmi = 0;
    for a = 1:size(joint_histogram, 1)
        for j = 1:size(joint_histogram, 2)
            qmi = qmi + (joint_histogram(a, j) - histogram_image1(a) * histogram_shifted_image2(j))^2;
        end
    end
    qmi_values(i) = qmi;
end

% Create a figure with two subplots to display results
figure;

% Plot 3: Correlation Coefficient vs. Shift (Negative of Image 1)
plot(shift_values, correlation_coefficients, '-o');
xlabel('Shift (tx pixels)');
ylabel('Correlation Coefficient');
title('Correlation Coefficient vs. Shift');

figure;
% Plot 4: QMI vs. Shift (Negative of Image 1)
plot(shift_values, qmi_values, '-o');
xlabel('Shift (tx pixels)');
ylabel('Quadratic Mutual Information (QMI)');
title('QMI vs. Shift');
