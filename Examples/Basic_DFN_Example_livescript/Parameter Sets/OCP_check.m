% To hard code a parameter set use the following line

% clear;
% load('negative_ocp.mat')
% 
% p = BPX_to_MATLAB('nmc_pouch_cell_BPX.json');
% 
% x_toofab = linspace(0, 1, 500);
% y_toofab = p.U_neg(x_toofab);
% 
% % Plot the function
% plot(x_toofab, y_toofab, 'r--', 'LineWidth', 1.5, 'DisplayName', 'TOOFAB');
% hold on
% plot(x_values, y_values, 'b-.', 'LineWidth', 1, 'DisplayName', 'PyBaMM');
% xlabel('State of Charge (SOC)');
% ylabel('OCP [V]');
% title('Negative Electrode Open Circuit Potential');
% hold off
% legend('Location', 'best'); 
% grid on;

clear;
load('negative_ocp.mat')

p = BPX_to_MATLAB('nmc_pouch_cell_BPX.json');

% Initialize variables
x_toofab = linspace(0, 1, 500); % Define the stoichiometric range
T = 298.15; % Temperature in Kelvin
y_toofab = zeros(size(x_toofab)); % Preallocate y_toofab for efficiency

% Loop to compute y_toofab values
for i = 1:length(x_toofab)
    y_toofab(i) = p.Ds_neg(x_toofab(i), T); % Compute negative electrode diffusivity
end

% Plot the function
plot(x_toofab, y_toofab, 'r--', 'LineWidth', 1.5, 'DisplayName', 'TOOFAB');
hold on
% plot(x_values, y_values, 'b-.', 'LineWidth', 1, 'DisplayName', 'PyBaMM');
xlabel('Stoich');
ylabel('Neg Diffusivity');
title('Negative Electrode Open Circuit Potential');
hold off
% legend('Location', 'best'); 
grid on;