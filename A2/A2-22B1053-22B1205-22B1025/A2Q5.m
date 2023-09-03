% Define the values and their corresponding probabilities
values = [1, 2, 3, 4, 5];
probabilities = [0.05, 0.4, 0.15, 0.3, 0.1];

% Values of N to consider
Ni = [5, 10, 20, 50, 100, 200, 500, 1000, 5000, 10000];

% Number of simulations for each N
nsamp = 6000;

% Initialize arrays to store MAD values and sample sizes
madValues = zeros(1, length(Ni));

% Loop through different values of N
for i = 1:length(Ni)
    N = Ni(i);
    
    % Simulate random variables for this N with the specified probabilities
    X = randsrc(N, nsamp, [values; probabilities]);
    
    % Calculate the average for each simulation
    X_avg = mean(X, 1);
    
    % Calculate the empirical CDF of X(N)
    [f, x] = ecdf(X_avg);
    
    % Calculate the mean and standard deviation (sigma) from X_avg
    mu = mean(X_avg);
    sigma = std(X_avg);
    
    % Calculate the CDF of a Gaussian with the same mean and variance
    f_gaussian = normcdf(x, mu, sigma);
    
    % Calculate the MAD between the empirical and Gaussian CDFs
    mad = max(abs(f - f_gaussian));
    
    % Store MAD value for this N
    madValues(i) = mad;
    
    % Print MAD value for this N
    fprintf("The MAD for N = %d is %f.\n", N, mad);
    
    % Create a new figure for Histogram
    figure;
    histogram(X_avg, 50, 'Normalization', 'probability');
    title(['Histogram of the Average, N = ', num2str(N)]);
    xlabel('Average Value');
    ylabel('Probability');
    
    % Save the histogram figure
    saveas(gcf, sprintf('N%d_histogram.png', N));
    
    % Create a new figure for ECDF and Gaussian CDF comparison
    figure;
    plot(x, f, '-b', 'LineWidth', 0.5);
    hold on;
    plot(x, f_gaussian, '-r', 'LineWidth', 0.5);
    xlabel(sprintf('Value, MAD = %f', mad));
    ylabel('Probability');
    legend('Empirical CDF', 'Gaussian CDF');
    title(['Empirical CDF vs. Gaussian CDF, N = ', num2str(N)]);
    
    % Save the CDF figure
    saveas(gcf, sprintf('N%d_cdf.png', N));
    
    % Pause to allow time for viewing the figures
    pause(5);
end

% Create a figure for MAD vs. N
figure;
plot(Ni, madValues, '-o', 'LineWidth', 2);
title('Mean Absolute Deviation (MAD) vs. Sample Size (N)');
xlabel('Sample Size (N)');
ylabel('MAD');
saveas(gcf, 'MAD_N.png');
