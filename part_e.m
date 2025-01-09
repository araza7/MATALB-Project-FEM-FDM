clear, clc;

% Parameters for FD
%Student ID 2024161152
A1 = 2;
A2 = 5;
L1 = 20 + A1 / 0.5; % Length of the first medium in meters
L2 = 15 + A2 / 0.5; % Length of the second medium in meters
c1 = 1500.0 + 10.0 * A1 - A2; % Speed of sound in the first medium in m/s
c2 = (1800.0 + 20.0 * A2 * A1) * (1 + 0.02 * 1i * (2 + sqrt(A1))); % Speed of sound in the second medium in m/s
rho = 1000; % Density in kg/m^3
f = 500.0 + A2 * 20.0 + A1; % Excitation frequency in Hz
w = 2 * pi * f; % Angular frequency

% Total length and number of points
L = L1 + L2;
num_discretizations = 10;
element_sizes = linspace(1.0, 0.1, num_discretizations); % From 1.0 m to 0.1 m
epsilon_FEM = zeros(num_discretizations, 1);
epsilon_FD = zeros(num_discretizations, 1);
npoints = zeros(num_discretizations, 1);

% Loop for FEM and FDM convergence analysis
for n = 1:num_discretizations
    dx = element_sizes(n); % Current element size
    N = round(L / dx); % Number of elements
    npoints(n) = N; % Store number of points for plotting

    % FEM Solution and Error Calculation
    P_FEM = fem_wave_solver(N - 1, L1, L2, c1, c2, rho, f); % FEM solver
    P_analytical = solution_1D(f, c1, c2, L1, L2, linspace(0, L, N)); % Analytical solution
    epsilon_FEM(n) = (sum(abs(P_analytical - P_FEM))) / N; % Calculate average absolute error for FEM

    % FDM Solution and Error Calculation
    P_FD = Helmholtz_Finite_Difference(N, L1, L2, c1, c2, f); % FDM solver
    epsilon_FD(n) = (sum(abs(P_analytical - P_FD))) / N; % Calculate average absolute error for FDM
end

% Plotting both convergence curves

figure;
hold on;
loglog(npoints+1, epsilon_FEM, 'o-', 'LineWidth', 2, 'DisplayName', 'FEM');
loglog(npoints+1,epsilon_FD, 's-', 'LineWidth', 2, 'DisplayName', 'FDM');
grid on;
xlabel('log_{10}(Number of Points, N)');
ylabel('log_{10}(Average Absolute Error, \epsilon)');
title('Convergence Analysis of Numerical Solutions');
legend('show');
set(gca,'xscale','log')
set(gca,'yscale','log')
hold off;

