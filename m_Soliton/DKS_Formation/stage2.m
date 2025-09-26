% --- MATLAB Code to Generate a Kerr Subcomb Schematic Diagram ---
% Author: Lu
% Date: 2025-09-25

% --- 1. Parameter Setup ---
clear; clc; close all;

FSR = 50;            % 
xi_modes = 11;        % Mode number between sub-combs and the primary comb (ξ)
delta_freq = 36;     % Spacing of sub-combs (δ)
num_cavity_modes = 30; % The number of cold-cavity resonance mode marks

% Powers
pump_power = 1.0;   % Normalized pump power
primary_power = 0.8;% Peak power for the primary comb
secondary_power = 0.4; % ~ 
tertiary_power = 0;  % ~

% Line number of sub-combs
subcomb_lines = 7; % ~
subcomb_mu = 2 * delta_freq; % The \mu-th spectral line to display

% --- 2. Initialization ---
figure('Color', 'w', 'Position', [100, 100, 1900, 1200]);% [100, 100, 1900, 700]
hold on;
ax = gca; 
ax.FontSize = 38;
ax.LineWidth = 1.2;
ax.TickDir = 'out';


% --- 3. Add the resonance frequencies of the cold cavity (Background) ---
cavity_mode_freqs = (-num_cavity_modes:num_cavity_modes) * FSR;
cavity_mode_power = 0.35 * ones(size(cavity_mode_freqs)); % the height of the cold cavity resonance mode
stem(cavity_mode_freqs, cavity_mode_power, ...
    'Color', [0.85, 0.85, 0.85], ... % light gray
    'Marker', 'none', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Cavity Modes');


% --- 4. Comb lines of sub-combs ---

% Define a function as the envelop of sub-combs
subcomb_envelope = @(x, x0, A, w) A * (w^2 ./ ((x - x0).^2 + w^2)); % Lorentzian envelop

% Pump (ξ=0)
pump_freqs = (-(subcomb_lines-1)/2:(subcomb_lines-1)/2) * delta_freq;
pump_powers = subcomb_envelope(pump_freqs, 0, pump_power, subcomb_mu);
stem(pump_freqs, pump_powers, ...
    'Color', 'r', 'Marker', 'none', 'LineWidth', 4, 'DisplayName', 'Pump');

% The 1st pair of sub-combs (±ξ)
primary_centers = [-1, 1] * xi_modes * FSR;
for center = primary_centers
    sub_freqs = center + (-(subcomb_lines-1)/2:(subcomb_lines-1)/2) * delta_freq;
    sub_powers = subcomb_envelope(sub_freqs, center, primary_power, subcomb_mu);
    stem(sub_freqs, sub_powers, ...
        'Color', [0, 0, 0.8], 'Marker', 'none', 'LineWidth', 4, 'DisplayName', 'Primary Subcombs');
end

% The 2nd pair of sub-combs (±2ξ)
secondary_centers = [-2, 2] * xi_modes * FSR;
for center = secondary_centers
    sub_freqs = center + (-(subcomb_lines-1)/2:(subcomb_lines-1)/2) * delta_freq;
    sub_powers = subcomb_envelope(sub_freqs, center, secondary_power, subcomb_mu);
    stem(sub_freqs, sub_powers, ...
        'Color', [0.4, 0.6, 1], 'Marker', 'none', 'LineWidth', 4, 'DisplayName', 'Secondary Subcombs');
end

% Higher order sub-combs
tertiary_centers = [-3, 3] * xi_modes * FSR;
for center = tertiary_centers
    sub_freqs = center + (-(subcomb_lines-1)/2:(subcomb_lines-1)/2) * delta_freq;
    sub_powers = subcomb_envelope(sub_freqs, center, tertiary_power, subcomb_mu);
    stem(sub_freqs, sub_powers, ...
        'Color', [0.7, 0.8, 1], 'Marker', 'none', 'LineWidth', 1.5);
end

% --- 5. Annotations ---
% Labels
xlabel('Frequency (\omega, a.u.)', 'FontSize', 44, 'FontWeight', 'bold');
ylabel('Power (P, a.u.)', 'FontSize', 44, 'FontWeight', 'bold');
% t = title('Stage 2: Sub-combs', 'FontSize', 24);
% 
% % Adjustment of title position
% current_pos = t.Position;
% t.Position = [current_pos(1), current_pos(2) + 0.22, current_pos(3)];


% The span of the axis
xlim([min(cavity_mode_freqs)-FSR, max(cavity_mode_freqs)+FSR])
ylim([0, 1.2]);
box on;

% Texts (ξ=0, ±ξ, etc.)
text(0, pump_power + 0.1, '\xi_0 = 0', 'FontSize', 36, 'HorizontalAlignment', 'center');
text(primary_centers(2), primary_power + 0.1, '\xi', 'FontSize', 36, 'HorizontalAlignment', 'center');
text(primary_centers(1), primary_power + 0.1, '-\xi', 'FontSize', 36, 'HorizontalAlignment', 'center');
text(secondary_centers(2), secondary_power + 0.15, '2\xi', 'FontSize', 36, 'HorizontalAlignment', 'center');
text(secondary_centers(1), secondary_power + 0.15, '-2\xi', 'FontSize', 36, 'HorizontalAlignment', 'center');


hold off;
legend('off'); % Exit current figure
