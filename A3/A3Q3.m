% Load XYZ coordinates from 'XYZ.txt'
data = dlmread('XYZ.txt');

% Separate X, Y, and Z coordinates
X = data(:, 1);
Y = data(:, 2);
Z = data(:, 3);

% Calculate the necessary summations
ZX = sum(X .* Z);
ZY = sum(Y .* Z);
X2 = sum(X .^ 2);
Y2 = sum(Y .^ 2);
sum_X = sum(X);
sum_Y = sum(Y);
sum_Z = sum(Z);
XY = sum(X .* Y);

% Create the coefficient matrix A of the linear equations
A = [X2, XY, sum_X; XY, Y2, sum_Y; sum_X, sum_Y, numel(X)];

% Create the column vector B of the right-hand side
B = [ZX; ZY; sum_Z];

% Solve the linear system A*x = B for the unknowns x
x = A \ B;

% Extract coefficients for the plane equation
a = x(1);
b = x(2);
c = x(3);

% Calculate the predicted equation of the plane
plane_equation = sprintf('Z = %.4f*X + %.4f*Y + %.4f', a, b, c);

% Calculate the predicted Z values for the plane
predicted_Z = a * X + b * Y + c;

% Calculate the residuals (noise)
residuals = Z - predicted_Z;

% Calculate the manually calculated noise variance
mean_residuals = sum(residuals) / numel(residuals);
squared_diff = (residuals - mean_residuals).^2;
noise_variance_manual = sum(squared_diff) / (numel(residuals) - 1);

% Display the predicted plane equation and manually calculated noise variance
fprintf('Predicted Plane Equation: %s\n', plane_equation);
fprintf('Manually Calculated Noise Variance: %.4f\n', noise_variance_manual);

% Plot the original data and the predicted plane
figure;
scatter3(X, Y, Z, 'b.');
hold on;
[xgrid, ygrid] = meshgrid(min(X):0.1:max(X), min(Y):0.1:max(Y));
zgrid = a * xgrid + b * ygrid + c;
mesh(xgrid, ygrid, zgrid, 'FaceAlpha', 0.5);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Fitted Plane to XYZ Data');
legend('Data Points', 'Fitted Plane');
grid on;
hold off;