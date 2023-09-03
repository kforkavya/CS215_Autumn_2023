x = [-3: 0.02: 3]; %Values from -3 to 3 at steps of 0.02
y = 6.5*sin(2.1*x + pi/3); %Applying function to my input x

n = numel(x); %Total no. of elements
f=0.6;
k = round(f*n); %Number of points to be experimented
indices_to_corrupt = randperm(n, k); %These are indices that will be corrupted
corrupt_values = 100 + (120 - 100) * rand(size(indices_to_corrupt)); %These are the values of the noises
z=y; %Copying y to z
z(indices_to_corrupt) = z(indices_to_corrupt) + corrupt_values; %Noise added to the graph

ymedian = zeros(size(z));
yquartile = zeros(size(z));
ymean = zeros(size(z));

for i = 1:numel(z)
    % Define the neighbourhood indices
    left_index = max(1, i - 8);
    right_index = min(numel(z), i + 8);
    neighbourhood = z(left_index:right_index);
    
    % Calculate the median of the neighbourhood
    ymedian(i) = median(neighbourhood);
    ymean(i) = mean(neighbourhood);
    yquartile(i) = prctile(neighbourhood, 25);
end

figure;
plot(x, y, 'LineWidth', 1, 'DisplayName', 'Original Sine Wave');
hold on;
plot(x, z, 'LineWidth', 1, 'DisplayName', 'Corrupted Sine Wave');
plot(x, ymedian, 'LineWidth', 1, 'DisplayName', 'Moving Median Filtered Sine Wave');
plot(x, ymean, 'LineWidth', 1, 'DisplayName', 'Moving Mean Filtered Sine Wave');
plot(x, yquartile, 'LineWidth', 1, 'DisplayName', 'Moving Quartile Filtered Sine Wave');
hold off;

title('Original, Corrupted, and Filtered Sine Waves');
xlabel('x');
ylabel('y');
legend('Location', 'best');
grid on;

% Calculating Error for moving median filtering
error_median = sum((y - ymedian).^2) / sum(y.^2);

% Calculating Error for moving mean filtering
error_mean = sum((y - ymean).^2) / sum(y.^2);

% Calculating Error for moving quartile filtering
error_quartile = sum((y - yquartile).^2) / sum(y.^2);

% Print Error values
fprintf('Error for Moving Median Filtering: %.4f\n', error_median);
fprintf('Error for Moving Mean Filtering: %.4f\n', error_mean);
fprintf('Error for Moving Quartile Filtering: %.4f\n', error_quartile);
