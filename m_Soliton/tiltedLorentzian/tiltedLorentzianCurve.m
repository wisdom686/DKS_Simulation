% --- MATLAB Code to Generate a Kerr Subcomb Schematic Diagram ---
% Author: Lu
% Date: 2025-09-26
clear; close all; clc;

% --- parameters ---
% f^2 : Normalized pump power, the larger it is, the more tilted the curve
S = 2;                

% Set the detuning range
Delta_vec = linspace(-6, 8, 2001);%6,8

% --- Solving cubic equation: I^3 - 2*Delta*I^2 + (1 + Delta^2)*I - S^2 = 0 ---
% roots_all: A num_points*3 matrix to store the roots of the function at sampled Delta points
num_points = numel(Delta_vec);
roots_all = NaN(num_points, 3); % store the three branches of the roots

% Calculate the roots at all detuning points
for k = 1:num_points
    D = Delta_vec(k);

    % coefficients of the cubic equation
    coeff = [1, -2*D, (1+D^2), -S^2];

    % solving the equation
    r = roots(coeff);

    % rule out negative roots
    r_real = r(abs(imag(r)) < 1e-8); % check if the imaginary part is close to 0
    r_pos = real(r_real);
    r_pos = r_pos(r_pos >= 0); % keep only non-negative roots

    % sort the roots
    r_pos = sort(r_pos, 'ascend');

    if k == 1
        % initialize at the first point
        roots_all(k, 1:numel(r_pos)) = r_pos;
    else
        prev_roots = roots_all(k-1, :);
        current_roots = NaN(1, 3);

        % Create mutable copies of roots to work with
        r_pos_available = r_pos; 
        prev_indices_unmatched = find(~isnan(prev_roots));

        % --- Globally Optimal Matching Algorithm ---
        % This loop iterates as many times as there are potential matches.
        % In each iteration, it finds the SINGLE BEST MATCH among all
        % remaining possibilities and assigns it.

        num_potential_matches = min(numel(prev_indices_unmatched), numel(r_pos_available));

        for match_iter = 1:num_potential_matches

            min_global_dist = inf;
            best_prev_idx = -1;
            best_r_pos_idx_in_available = -1;

            % Iterate through all currently unmatched previous branches
            for p_idx = prev_indices_unmatched
                % Find the closest root in the currently available r_pos list
                [min_local_dist, local_idx] = min(abs(r_pos_available - prev_roots(p_idx)));

                % If this pair is better than the best found so far, record it
                if min_local_dist < min_global_dist
                    min_global_dist = min_local_dist;
                    best_prev_idx = p_idx;
                    best_r_pos_idx_in_available = local_idx;
                end
            end

            % After checking all possibilities, we have the globally best pair.
            % Now, make the assignment.
            if best_prev_idx ~= -1
                % Assign the value to the correct branch in current_roots
                current_roots(best_prev_idx) = r_pos_available(best_r_pos_idx_in_available);

                % Remove the assigned root and branch from the available lists
                r_pos_available(best_r_pos_idx_in_available) = [];
                prev_indices_unmatched(prev_indices_unmatched == best_prev_idx) = [];
            else
                % This should not happen if there are matches to be made, but as a safeguard:
                break;
            end
        end

        % --- Handle New Branches ---
        % If there are any roots left in r_pos_available, they are new branches.
        if ~isempty(r_pos_available)
            nan_slots = find(isnan(current_roots));
            num_to_fill = min(numel(r_pos_available), numel(nan_slots));
            % Assign them to the empty slots, sorted by value for consistency
            current_roots(nan_slots(1:num_to_fill)) = sort(r_pos_available(1:num_to_fill));
        end

        % Store the correctly tracked roots
        roots_all(k, :) = current_roots;
    end
end


% --- Plotting ---
figure('Position', [100, 100, 1900, 1200]);  
hold on;

% plot the linear Lorentzian curve as a reference
I_linear = S^2 ./ (1 + Delta_vec.^2);
plot(Delta_vec, I_linear, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'DisplayName', 'Linear Cavity');


% Find jump point indices before plotting (ensure they are calculated for truncation)
minD = min(Delta_vec);
maxD = max(Delta_vec);
bp1_idx = find(~isnan(roots_all(:,3)), 1, 'first');
bp2_idx = find(~isnan(roots_all(:,3)), 1, 'last');
bp1 = bp1_idx*(maxD-minD)/numel(Delta_vec)+minD;
bp2 = bp2_idx*(maxD-minD)/numel(Delta_vec)+minD;
% xline(bp1, 'm--', 'LineWidth', 2,'HandleVisibility','off');
% xline(bp2, 'm--', 'LineWidth', 2,'HandleVisibility','off');
h = line([bp1 bp1], [roots_all(bp1_idx,1)  roots_all(bp1_idx,2)], 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 1.5);
h.Annotation.LegendInformation.IconDisplayStyle = 'off'; % do not display in legend
h = line([bp2 bp2], [roots_all(bp2_idx,1)  roots_all(bp2_idx,2)], 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 1.5);
h.Annotation.LegendInformation.IconDisplayStyle = 'off'; % do not display in legend

% plot the three branches, using different styles to distinguish them
% stable lower branch (usually roots_all(:,1))
plot(Delta_vec, roots_all(:,1), 'b', 'LineWidth', 2, 'DisplayName', 'Upper Branch');      
% stable upper branch (usually roots_all(:,3))
plot(Delta_vec, roots_all(:,3), 'r:', 'LineWidth', 2, 'DisplayName', 'Unstable Branch');      
% middle (unstable) branch (usually roots_all(:,2))
plot(Delta_vec, roots_all(:,2), 'r', 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'DisplayName', 'Lower Branch');  

grid on; 
box on;
xlabel('$\zeta_0$ (a.u.)','Interpreter','latex'); 
ylabel('$|\Psi|^2$ (a.u.)','Interpreter','latex');
title(sprintf('Tilted Lorentzian Curve, $f^2$ = %.2f', S),'Interpreter','latex');
set(gca, 'FontSize', 38); % adjust font size for readability

xlim([min(Delta_vec) max(Delta_vec)]);
ylim([0, max(roots_all(:)) * 1.1]);
legend('Location', 'NorthWest','FontName','Times New Roman'); % use DisplayName for legend
hold off;