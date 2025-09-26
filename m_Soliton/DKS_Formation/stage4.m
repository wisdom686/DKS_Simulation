% --- MATLAB Code to Generate a DKS Formation Schematic ---
% --- Stage 4: Dissipative Kerr Soliton (sech^2 envelope) ---
% Author: Lu
% Date: 2025-09-25

% --- 1. Parameter Setup ---
clear; clc; close all;
FSR = 50;           % Free Spectral Range of the cold cavity modes
num_cavity_modes = 25; % Number of comb lines to plot (on each side of the pump)

% DKS Envelope Parameters
peak_power = 0.9;   % Peak power of the sech^2 envelope
sech_width_modes = 10; % The characteristic width of the envelope, in units of FSR
                       % A smaller value means a broader spectrum

% --- 2. Initialization ---
figure('Color', 'w', 'Position', [100, 100, 1900, 1200]);
hold on;
ax = gca; 
ax.FontSize = 38;
ax.LineWidth = 1.2;
ax.TickDir = 'out';

% --- 3. Add the resonance frequencies of the cold cavity (Background) ---
% For the final DKS stage, we often omit the background modes for clarity,
% but I've included the code here, commented out, in case you want it.
% cavity_mode_freqs_bg = (-num_cavity_modes-10:num_cavity_modes+10) * FSR;
% cavity_mode_power_bg = 0.35 * ones(size(cavity_mode_freqs_bg));
% stem(cavity_mode_freqs_bg, cavity_mode_power_bg, ...
%     'Color', [0.85, 0.85, 0.85], ...
%     'Marker', 'none', ...
%     'LineWidth', 1.5);

% --- 4. Comb lines of the Dissipative Kerr Soliton ---

% --- Step 4a: Generate the comb line frequencies ---
% In the stable soliton state, the comb lines are evenly spaced by FSR.
comb_freqs = (-num_cavity_modes:num_cavity_modes) * FSR;

% --- Step 4b: Define the sech^2 envelope function ---
% The spectral envelope of a soliton is proportional to sech^2.
sech_width_freqs = sech_width_modes * FSR; % Convert width to frequency units
sech2_envelope = @(x, A, width) A * (sech(x / width)).^2;

% --- Step 4c: Calculate the power of each comb line ---
comb_powers = sech2_envelope(comb_freqs, peak_power, sech_width_freqs);

% --- Step 4d: Plot the comb lines (the blue stems) ---
stem(comb_freqs, comb_powers, ...
    'Color', [0, 0, 0.7], 'Marker', 'none', 'LineWidth', 4);

% --- Step 4e: Plot the smooth red envelope line over the top ---
% Create a much finer frequency grid for a smooth plot
fine_freq_grid = linspace(min(comb_freqs), max(comb_freqs), 1000);
envelope_line_power = sech2_envelope(fine_freq_grid, peak_power, sech_width_freqs);
plot(fine_freq_grid, envelope_line_power, 'r-', 'LineWidth', 4);

% --- Step 4f (Optional): Make the center line slightly taller, as in the reference ---
center_idx = find(comb_freqs == 0);
center_line_power = sech2_envelope(0, peak_power, sech_width_freqs) * 1.1; % 10% taller
stem(0, center_line_power, 'Color', [0, 0, 0.7], 'Marker', 'none', 'LineWidth', 4);


% --- 5. Annotations ---
% Add the new text annotations
text(num_cavity_modes*FSR*0.4, peak_power*0.5, 'sech^2 envelop', ...
    'Color', 'r', 'FontSize', 44, 'FontWeight', 'bold');
text(num_cavity_modes*FSR*0.01, 0.95, '\leftarrow pump', ...
    'FontSize', 44, 'FontWeight', 'bold');

% Labels
xlabel('Frequency (\omega, a.u.)', 'FontSize', 44, 'FontWeight', 'bold');
ylabel('Power (P, a.u.)', 'FontSize', 44, 'FontWeight', 'bold');
% t = title('Stage 4: Dissipative Kerr Soliton', 'FontSize', 24);
% 
% % Adjustment of title position
% current_pos = t.Position;
% t.Position = [current_pos(1), current_pos(2) + 0.22, current_pos(3)];

% The span of the axis
xlim([-num_cavity_modes*FSR - 2*FSR, num_cavity_modes*FSR + 2*FSR]);
ylim([0, 1.2]);
box on;
hold off;
legend('off'); % Exit current figure



% % --- 5. Annotations ---
% % As requested, this section is minimal.
% % Labels
% xlabel('Frequency (\omega)', 'FontSize', 24, 'FontWeight', 'bold');
% ylabel('Power (P)', 'FontSize', 24, 'FontWeight', 'bold');
% t = title('Stage 4: Dissipative Kerr Soliton', 'FontSize', 24);
% 
% % Adjustment of title position
% current_pos = t.Position;
% t.Position = [current_pos(1), current_pos(2) + 0.22, current_pos(3)];
% 
% % The span of the axis
% xlim([-num_cavity_modes*FSR - 2*FSR, num_cavity_modes*FSR + 2*FSR]);
% ylim([0, 1.2]);
% box on;
% 
% hold off;
% legend('off'); % Exit current figure