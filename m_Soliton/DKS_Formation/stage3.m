% --- MATLAB Code to Generate a Kerr Comb 'Merging' Stage Schematic ---
% --- Stage 3: Merging ---
% Author: Lu
% Date: 2025-09-25

% --- 1. Parameter Setup ---
clear; clc; close all;
FSR = 50;           % Free Spectral Range of the cold cavity modes
num_cavity_modes = 30; % The number of cold-cavity resonance mode marks

% --- Parameters for the prominent peaks (from your Stage 2 script) ---
xi_modes = 11;        % Mode number between sub-combs and the primary comb (ξ)
pump_power = 1.0;     % Peak power for the pump line
primary_power = 0.8;  % Peak power for the primary sub-comb centers
secondary_power = 0.4;% Peak power for the secondary sub-comb centers

% --- Parameters for the merged/chaotic background spectrum ---
base_power = 0.1;     % The minimum power level for the background comb lines
max_wavy_power = 0.6; % The maximum power the wavy envelope can reach
noise_level = 0.2;    % The amount of randomness/noise

% --- 2. Initialization ---
figure('Color', 'w', 'Position', [100, 100, 1900, 1200]);
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
    'LineWidth', 1.5);

% --- 4. Comb lines of the Merging Stage ---

% --- Step 4a: Generate the low-power, wavy, chaotic background spectrum ---
all_comb_freqs = cavity_mode_freqs;
% Define a wavy envelope to modulate the power
wavy_envelope = @(x) base_power + (max_wavy_power - base_power) / 2 * (1 + cos(2 * pi * x / (xi_modes * FSR)));
all_comb_powers = wavy_envelope(all_comb_freqs);
% Add randomness to make it look more chaotic
random_modulation = 1 - noise_level * rand(size(all_comb_powers));
all_comb_powers = all_comb_powers .* random_modulation;

% --- Step 4b: Manually set the power of prominent peaks to Stage 2 levels ---
% Define the center frequencies of the prominent peaks
pump_center_freq = 0;
primary_center_freqs = [-1, 1] * xi_modes * FSR;
secondary_center_freqs = [-2, 2] * xi_modes * FSR;

% Find the indices of these frequencies in our frequency array
pump_idx = find(all_comb_freqs == pump_center_freq);
primary_idx1 = find(all_comb_freqs == primary_center_freqs(1));
primary_idx2 = find(all_comb_freqs == primary_center_freqs(2));
secondary_idx1 = find(all_comb_freqs == secondary_center_freqs(1));
secondary_idx2 = find(all_comb_freqs == secondary_center_freqs(2));

% Override the power at these specific indices
all_comb_powers(pump_idx) = pump_power;
all_comb_powers([primary_idx1, primary_idx2]) = primary_power;
all_comb_powers([secondary_idx1, secondary_idx2]) = secondary_power;

% --- Step 4c: Plot the comb lines with appropriate coloring ---
center_width_modes = 4; 
primary_width_modes = 15;

for i = 1:length(all_comb_freqs)
    freq = all_comb_freqs(i);
    power = all_comb_powers(i);
    mode_index = abs(round(freq / FSR));
    
    line_color = [0.5, 0.2, 0.8]; % Default color (purple)
    
    if mode_index <= center_width_modes
        line_color = 'r'; % Red for the center
    elseif mode_index <= primary_width_modes
        line_color = [0, 0, 0.8]; % Dark blue for the primary area
    end
    
    stem(freq, power, 'Color', line_color, 'Marker', 'none', 'LineWidth', 4);
end

% --- 5. Annotations ---
% Labels
xlabel('Frequency (\omega, a.u.)', 'FontSize', 44, 'FontWeight', 'bold');
ylabel('Power (P, a.u.)', 'FontSize', 44, 'FontWeight', 'bold');
% t = title('Stage 3: Merging', 'FontSize', 24);
% 
% % Adjustment of title position
% current_pos = t.Position;
% t.Position = [current_pos(1), current_pos(2) + 0.22, current_pos(3)];

% The span of the axis
xlim([min(cavity_mode_freqs)-FSR, max(cavity_mode_freqs)+FSR]);
ylim([0, 1.2]);
box on;

hold off;
legend('off'); % Exit current figure
