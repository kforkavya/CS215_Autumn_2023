clear;
clc;

rng(0);

% Generate n = 1000 independent samples from N(0, 16)
n = 1000;
mu = 0;
std_dev = sqrt(16); % Standard Deviation
data = mu + std_dev * randn(n, 1);

% Split the data into training set T (750 samples) ...
random_indices = randperm(n);
T_size = 750;
T_indices = random_indices(1:T_size);
% Create the training set T and validation set V
T = data(T_indices);

% Values of sigma to evaluate
sigma_values = [0.001, 0.1, 0.2, 0.9, 1, 2, 3, 5, 10, 20, 100];

% Initialize arrays to store LL values, D values and log(sigma) values
LL_values = zeros(1, numel(sigma_values));
D_values = zeros(1, numel(sigma_values));
log_sigma_values = log(sigma_values);

% Initialize the best_sigma and best_LL
best_sigma = 0;
best_LL = -inf;
best_D = inf;

% Calculate LL for each value of sigma and find the best_sigma
for i = 1:numel(sigma_values)
    sigma_bandwidth = sigma_values(i);
    
    % Calculate KDE for each data point in V using sigma
    likelihood = 0;
    D = 0;
    for T2_index = 1:T_size
        vi = T(T2_index);
        kernel = 0;
        for T_index = 1:T_size
            ti = T(T_index);
            kernel = kernel + (1 / (T_size * sigma_bandwidth * sqrt(2 * pi))) * exp(-(vi - ti)^2 / (2 * sigma_bandwidth^2));
        end
        true_density_vi = 1 / (std_dev * sqrt(2 * pi)) * exp(-(vi - mu).^2 / (2 * std_dev^2));
        D = D + (true_density_vi - kernel)^2;
        likelihood = likelihood + log(kernel);
    end
    
    LL_values(i) = likelihood;
    D_values(i) = D;
    
    % Check if this sigma gives a higher likelihood
    if likelihood > best_LL
        best_LL = likelihood;
        best_sigma = sigma_bandwidth;
        best_LL_D_value = D; %This is D value of best LL
    end

    % Checking for sigma with minimum D value
    if D < best_D
        best_D = D;
        sigma_best_D = sigma_bandwidth;
    end
end

% Plot LL versus log(sigma)
figure;
plot(log_sigma_values, LL_values, '-o');
xlabel('log(sigma)');
ylabel('Log Joint Likelihood (LL)');
title('Log Joint Likelihood vs. log(sigma) for T=V');
grid on;

% Plot D versus log(sigma)
figure;
plot(log_sigma_values, D_values, '-o');
xlabel('log(sigma)');
ylabel('D Value (D)');
title('D Value vs. log(sigma) for T=V');
grid on;

% Display the best sigma and LL value and D values
fprintf('Best Sigma for Max. LL for T=V: %.4f\n', best_sigma);
fprintf('Best Log Joint Likelihood (LL) for T=V: %.4f\n', best_LL);
fprintf('D Value for Best LL for T=V: %.4f\n\n', best_LL_D_value);
fprintf('Best Sigma for Min. D for T=V: %.4f\n', sigma_best_D);
fprintf('Best D Value for T=V: %.4f\n', best_D);

% Calculate p_hat values for the best sigma
x_values = -8:0.1:8;
p_hat_values = zeros(numel(x_values), 1);

% Estimate the density pˆn(x;σ) for x ∈ [-8 : 0.1 : 8] using the best_sigma
% with maximum log likelihood
for j = 1:numel(x_values)
    x = x_values(j);
    p_hat = 0;
    for T_index = 1:T_size
        ti = T(T_index);
        p_hat = p_hat + (1 / (T_size * best_sigma * sqrt(2 * pi))) * exp(-(x - ti)^2 / (2 * best_sigma^2));
    end
    p_hat_values(j) = p_hat;
end

% Plot estimated density for the best sigma and overlay true density
true_density = 1 / (std_dev * sqrt(2 * pi)) * exp(-(x_values - mu).^2 / (2 * std_dev^2));

figure;
plot(x_values, p_hat_values, 'b', 'LineWidth', 2, 'DisplayName', 'Estimated Density');
hold on;
plot(x_values, true_density, 'r', 'LineWidth', 2, 'DisplayName', 'True Density');
xlabel('x');
ylabel('Density');
title('Estimated Density vs. True Density (Best LL) for T=V');
legend('Location', 'North');
grid on;

p_hat_values = zeros(numel(x_values), 1);

% Estimate the density pˆn(x;σ) for x ∈ [-8 : 0.1 : 8] using the best_sigma
% with minimum D Value
for j = 1:numel(x_values)
    x = x_values(j);
    p_hat = 0;
    for T_index = 1:T_size
        ti = T(T_index);
        p_hat = p_hat + (1 / (T_size * sigma_best_D * sqrt(2 * pi))) * exp(-(x - ti)^2 / (2 * sigma_best_D^2));
    end
    p_hat_values(j) = p_hat;
end

figure;
plot(x_values, p_hat_values, 'b', 'LineWidth', 2, 'DisplayName', 'Estimated Density');
hold on;
plot(x_values, true_density, 'r', 'LineWidth', 2, 'DisplayName', 'True Density');
xlabel('x');
ylabel('Density');
title('Estimated Density vs. True Density (Best/Min D) for T=V');
legend('Location', 'North');
grid on;
