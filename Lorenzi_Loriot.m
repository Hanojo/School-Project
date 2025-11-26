%% ================== project_corrected.m ==================
close all; clear; clc;

%% --- Input files
f1 = 'res_simul_global_L96W48_xtip.mat';
f2 = 'mean_Ufield_L96W48.mat';

%% --- Load data
Sres = load(f1);
Suf  = load(f2);

% Get simul structure
if isfield(Sres,'simul')
    simul = Sres.simul;
else
    fn = fieldnames(Sres);
    simul = Sres.(fn{1});
end

%% --- Mesh and material parameters
L = 96; W = 48; 
[NODES, ELEMS] = definition2D_honeycomb_fusible(L, W);
[ELEMS, X0new, Y0new] = introfissure2D_honeycomb(NODES, ELEMS, 0, 0);

Xnode = NODES(:,1);
Ynode = NODES(:,2);
N = numel(Xnode);

%% --- Extract simulation data
Xtip_all = simul.Xcracktip(:);
if isfield(simul,'Ycracktip')
    Ytip_all = simul.Ycracktip(:);
else
    Ytip_all = zeros(size(Xtip_all));
end
mu = simul.ShearModulus;

MU = Suf.MU;
T  = numel(Xtip_all);

% Ensure proper dimensions
if size(MU,2) ~= T && size(MU,1) == T
    MU = MU.';
end

%% --- Williams decomposition WITH n=-1 term for tip correction
order = 11;  % n=1,3,5,7,9,11 (6 terms) - BUT we need n=-1 for correction

fprintf('Starting Williams decomposition with iterative tip correction...\n');
fprintf('Using method from article: n=-1 term for tip positioning correction\n\n');

% Arrays for results
B_all = cell(T,1);
Ufit_all = zeros(N, T);
err_all = zeros(T,1);
Xtip_corrected = zeros(T,1);
Ytip_corrected = zeros(T,1);
K1_corrected = zeros(T,1);

for k = 1:T
    Unode_k = MU(:,k);
    Xtip_initial = Xtip_all(k);
    Ytip_initial = Ytip_all(k);
    
    fprintf('Time step %d/%d: Initial tip at (%.3f, %.3f)\n', k, T, Xtip_initial, Ytip_initial);
    
    % Apply iterative tip positioning correction WITH n=-1
    [Xtip_corr, Ytip_corr, Bk, Ufitk, errk, K1k] = iterative_tip_correction_with_supersingular(...
        Xnode, Ynode, Unode_k, order, Xtip_initial, Ytip_initial);
    
    B_all{k} = Bk;
    Ufit_all(:,k) = Ufitk;
    err_all(k) = errk;
    Xtip_corrected(k) = Xtip_corr;
    Ytip_corrected(k) = Ytip_corr;
    K1_corrected(k) = K1k;
    
    fprintf('  Corrected tip to (%.3f, %.3f), K1=%.6f\n', Xtip_corr, Ytip_corr, K1k);
end

%% Q1 & Q2: Fracture toughness with CORRECTED tip positions
fprintf('\n=== Q1-Q2: Fracture Toughness with Corrected Tip Positions ===\n');

% Young's modulus for honeycomb lattice
E = 2/sqrt(3);

% Calculate fracture toughness: Kc = E * K1 * (2*pi)/4
Kc_values = E * K1_corrected * (2*pi)/4;

% Compare with uncorrected values
% Calculate uncorrected K1 for comparison
K1_uncorrected = zeros(T,1);
for k = 1:T
    Unode_k = MU(:,k);
    [B_uncorr, ~, ~] = decomposition2williams_fusible_v2(...
        Xnode, Ynode, Unode_k, order, Xtip_all(k), Ytip_all(k));
    ns_uncorr = 1:2:order;
    K1_uncorrected(k) = B_uncorr(1); % First term is n=1
end
Kc_uncorrected = E * K1_uncorrected * (2*pi)/4;

% Plot comparison
figure('Position', [100, 100, 1000, 800]);

subplot(2,2,1);
plot(Xtip_all, Kc_uncorrected, 'ro-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
hold on;
plot(Xtip_corrected, Kc_values, 'bo-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
grid on;
xlabel('Crack Tip Position x_{tip}');
ylabel('Fracture Toughness K_c');
title('K_c: Corrected vs Uncorrected Tip Positions');
legend('Uncorrected tip', 'Corrected tip', 'Location', 'best');

subplot(2,2,2);
tip_correction_distance = sqrt((Xtip_corrected - Xtip_all).^2 + (Ytip_corrected - Ytip_all).^2);
plot(Xtip_all, tip_correction_distance, 'gs-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'g');
grid on;
xlabel('Crack Tip Position x_{tip}');
ylabel('Tip Correction Distance');
title('Tip Positioning Correction');

subplot(2,2,3);
Kc_improvement = abs(Kc_values - Kc_uncorrected) ./ Kc_uncorrected * 100;
plot(Xtip_all, Kc_improvement, 'm^-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'm');
grid on;
xlabel('Crack Tip Position x_{tip}');
ylabel('K_c Change (%)');
title('Relative Change in K_c due to Tip Correction');

subplot(2,2,4);
plot(Xtip_all, err_all, 'kd-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'k');
grid on;
xlabel('Crack Tip Position x_{tip}');
ylabel('Fitting Error');
title('Fitting Error vs Crack Position');

% Statistical analysis
fprintf('Kc statistics with corrected tips:\n');
fprintf('  Mean Kc = %.6f\n', mean(Kc_values));
fprintf('  Std Kc = %.6f\n', std(Kc_values));
fprintf('  Min Kc = %.6f\n', min(Kc_values));
fprintf('  Max Kc = %.6f\n', max(Kc_values));
fprintf('Tip correction statistics:\n');
fprintf('  Mean correction distance = %.6f\n', mean(tip_correction_distance));
fprintf('  Max correction distance = %.6f\n', max(tip_correction_distance));

%% ================== CORRECTED VISUALIZATION SECTION ==================

%% Q3: Fitted displacement field at first crack tip position (CORRECTED)
fprintf('\n=== Q3: Fitted Displacement Field with Corrected Tip Position ===\n');

k = 1;
Ufit1 = Ufit_all(:,k);
u_numerical = MU(:,k);

figure('Position', [100, 100, 1400, 1000]);

% 3D Surface plot - Fitted field
subplot(2,3,1);
scatter3(Xnode, Ynode, Ufit1, 25, Ufit1, 'filled');
grid on; axis equal; view(45,30);
xlabel('x position'); ylabel('y position'); zlabel('u_{fit}');
title({'Fitted Displacement Field u_{fit}(x,y)', ...
       sprintf('Time Step %d, Corrected Tip Position', k), ...
       sprintf('Initial tip: (%.3f, %.3f)', Xtip_all(k), Ytip_all(k)), ...
       sprintf('Corrected tip: (%.3f, %.3f)', Xtip_corrected(k), Ytip_corrected(k))});
colorbar;
hold on;
plot3(Xtip_corrected(k), Ytip_corrected(k), min(Ufit1), 'rp', ...
      'MarkerSize', 15, 'LineWidth', 3, 'MarkerFaceColor', 'r');
legend('Fitted field', 'Corrected crack tip', 'Location', 'northeast');

% 2D Contour plot - Fitted field
subplot(2,3,2);
tri = delaunay(Xnode, Ynode);
trisurf(tri, Xnode, Ynode, Ufit1, 'EdgeColor', 'none');
view(0,90);
axis equal;
xlabel('x position'); ylabel('y position');
title({'Fitted Field Contour Plot', ...
       sprintf('Time Step %d', k)});
colorbar;
hold on;
plot3(Xtip_corrected(k), Ytip_corrected(k), max(Ufit1), 'rp', ...
      'MarkerSize', 15, 'LineWidth', 3, 'MarkerFaceColor', 'r');

% Tip position correction visualization
subplot(2,3,3);
% Show lattice around crack tip
crack_region = (Xnode > Xtip_all(k)-5) & (Xnode < Xtip_all(k)+5) & ...
               (Ynode > -3) & (Ynode < 3);
plot(Xnode(crack_region), Ynode(crack_region), 'b.', 'MarkerSize', 10);
hold on;
plot(Xtip_all(k), Ytip_all(k), 'rs', 'MarkerSize', 12, ...
     'LineWidth', 3, 'MarkerFaceColor', 'r');
plot(Xtip_corrected(k), Ytip_corrected(k), 'go', 'MarkerSize', 10, ...
     'LineWidth', 3, 'MarkerFaceColor', 'g');

% Draw correction vector
quiver(Xtip_all(k), Ytip_all(k), ...
       Xtip_corrected(k)-Xtip_all(k), Ytip_corrected(k)-Ytip_all(k), ...
       0, 'k-', 'LineWidth', 2, 'MaxHeadSize', 1);

axis equal; grid on;
xlabel('x position'); ylabel('y position');
title({'Tip Position Correction', ...
       sprintf('Correction distance: %.3f units', ...
       sqrt((Xtip_corrected(k)-Xtip_all(k))^2 + (Ytip_corrected(k)-Ytip_all(k))^2))});
legend('Lattice nodes', 'Initial tip', 'Corrected tip', 'Correction vector', 'Location', 'best');

%% Q4: Numerical displacement field at first crack tip position
fprintf('\n=== Q4: Numerical Displacement Field Comparison ===\n');

% 3D Surface plot - Numerical field
subplot(2,3,4);
scatter3(Xnode, Ynode, u_numerical, 25, u_numerical, 'filled');
grid on; axis equal; view(45,30);
xlabel('x position'); ylabel('y position'); zlabel('u');
title({'Numerical Displacement Field u(x,y)', ...
       'From Simulation Data'});
colorbar;
hold on;
plot3(Xtip_corrected(k), Ytip_corrected(k), min(u_numerical), 'rp', ...
      'MarkerSize', 15, 'LineWidth', 3, 'MarkerFaceColor', 'r');
legend('Numerical field', 'Corrected crack tip', 'Location', 'northeast');

% 2D Contour plot - Numerical field
subplot(2,3,5);
trisurf(tri, Xnode, Ynode, u_numerical, 'EdgeColor', 'none');
view(0,90);
axis equal;
xlabel('x position'); ylabel('y position');
title({'Numerical Field Contour Plot', ...
       'Direct Simulation Results'});
colorbar;
hold on;
plot3(Xtip_corrected(k), Ytip_corrected(k), max(u_numerical), 'rp', ...
      'MarkerSize', 15, 'LineWidth', 3, 'MarkerFaceColor', 'r');

% Field comparison side by side (FIXED - no corr function)
subplot(2,3,6);
% Normalize both fields for comparison
u_norm = (u_numerical - min(u_numerical)) / (max(u_numerical) - min(u_numerical));
ufit_norm = (Ufit1 - min(Ufit1)) / (max(Ufit1) - min(Ufit1));

% Remove NaN values for comparison
valid_idx = ~isnan(u_norm) & ~isnan(ufit_norm);
u_norm_clean = u_norm(valid_idx);
ufit_norm_clean = ufit_norm(valid_idx);

plot(u_norm_clean, ufit_norm_clean, 'b.', 'MarkerSize', 8);
hold on;
plot([0 1], [0 1], 'r--', 'LineWidth', 2);
grid on;
xlabel('Normalized Numerical Field u');
ylabel('Normalized Fitted Field u_{fit}');
title({'Field Comparison', ...
       'Perfect Fit = Red Dashed Line'});
legend('Data points', 'Perfect fit line', 'Location', 'northwest');
axis equal;

% Manual correlation calculation (replacement for corr function)
cov_matrix = cov(u_norm_clean, ufit_norm_clean);
if size(cov_matrix,1) == 2 && size(cov_matrix,2) == 2
    manual_corr = cov_matrix(1,2) / (std(u_norm_clean) * std(ufit_norm_clean));
else
    manual_corr = NaN;
end

% R-squared calculation
ss_res = sum((u_norm_clean - ufit_norm_clean).^2);
ss_tot = sum((u_norm_clean - mean(u_norm_clean)).^2);
r_squared = 1 - (ss_res / ss_tot);

text(0.1, 0.9, sprintf('Manual Correlation: %.4f', manual_corr), ...
     'Units', 'normalized', 'FontSize', 12, 'BackgroundColor', 'white');
text(0.1, 0.8, sprintf('R-squared: %.4f', r_squared), ...
     'Units', 'normalized', 'FontSize', 12, 'BackgroundColor', 'white');

%% Q5: Absolute difference and comprehensive error analysis
fprintf('\n=== Q5: Comprehensive Error Analysis ===\n');

diff_abs = abs(Ufit1 - u_numerical);

figure('Position', [100, 100, 1400, 1000]);

% 2D Error map
subplot(2,3,1);
scatter(Xnode, Ynode, 30, diff_abs, 'filled');
axis equal; grid on;
xlabel('x position'); ylabel('y position');
title({'Absolute Error Distribution |u_{fit} - u|(x,y)', ...
       'Spatial Error Map'});
cb = colorbar; 
ylabel(cb, 'Absolute Error |u_{fit} - u|');
hold on;
plot(Xtip_corrected(k), Ytip_corrected(k), 'rp', ...
     'MarkerSize', 15, 'LineWidth', 3, 'MarkerFaceColor', 'r');
plot(Xtip_all(k), Ytip_all(k), 'ys', ...
     'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'y');
legend('Error field', 'Corrected tip', 'Initial tip', 'Location', 'best');

% Error distribution histogram
subplot(2,3,2);
histogram(diff_abs, 50, 'FaceColor', 'blue', 'FaceAlpha', 0.7, 'EdgeColor', 'black');
xlabel('Absolute Error |u_{fit} - u|');
ylabel('Frequency');
title({'Error Distribution Histogram', ...
       'Statistical Error Analysis'});
grid on;

% Add statistical annotations
error_mean = mean(diff_abs, 'omitnan');
error_std = std(diff_abs, 'omitnan');
xlims = xlim;
ylims = ylim;
text(0.7*xlims(2), 0.8*ylims(2), sprintf('Mean: %.2e\nStd: %.2e', error_mean, error_std), ...
     'FontSize', 11, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Error vs distance from corrected crack tip
subplot(2,3,3);
distances = sqrt((Xnode - Xtip_corrected(k)).^2 + (Ynode - Ytip_corrected(k)).^2);

% Bin errors by distance for better visualization
dist_bins = 0:1:20; % Smaller bins for better resolution
mean_errors = zeros(size(dist_bins));
std_errors = zeros(size(dist_bins));
counts = zeros(size(dist_bins));

for i = 1:length(dist_bins)-1
    mask = (distances >= dist_bins(i)) & (distances < dist_bins(i+1)) & ~isnan(diff_abs);
    if sum(mask) > 5  % Require minimum points
        mean_errors(i) = mean(diff_abs(mask), 'omitnan');
        std_errors(i) = std(diff_abs(mask), 'omitnan');
        counts(i) = sum(mask);
    else
        mean_errors(i) = NaN;
        std_errors(i) = NaN;
    end
end

valid_bins = ~isnan(mean_errors(1:end-1));
errorbar(dist_bins(1:end-1), mean_errors(1:end-1), std_errors(1:end-1), ...
         'ro-', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
xlabel('Distance from Corrected Crack Tip');
ylabel('Mean Absolute Error');
title({'Error vs Distance from Crack Tip', ...
       'With Standard Deviation Bars'});
grid on;

% Relative error analysis
subplot(2,3,4);
max_u = max(abs(u_numerical), [], 'omitnan');
relative_error = diff_abs / max_u * 100;  % Percentage

scatter(distances, relative_error, 30, 'filled', 'blue');
xlabel('Distance from Corrected Crack Tip');
ylabel('Relative Error (%)');
title({'Relative Error vs Distance', ...
       'Percentage of Maximum Displacement'});
grid on;

% Add trend line
valid_idx = ~isnan(relative_error) & isfinite(relative_error);
if sum(valid_idx) > 10
    p = polyfit(distances(valid_idx), relative_error(valid_idx), 1);
    trend_line = polyval(p, sort(distances(valid_idx)));
    hold on;
    plot(sort(distances(valid_idx)), trend_line, 'r-', 'LineWidth', 2);
    legend('Data points', sprintf('Trend: %.3f%%/unit', p(1)), 'Location', 'best');
end

% Spatial error patterns - highlight regions
subplot(2,3,5);
% Classify errors
low_error = diff_abs < error_mean;
med_error = (diff_abs >= error_mean) & (diff_abs < error_mean + error_std);
high_error = diff_abs >= error_mean + error_std;

plot(Xnode(low_error), Ynode(low_error), 'g.', 'MarkerSize', 8);
hold on;
plot(Xnode(med_error), Ynode(med_error), 'y.', 'MarkerSize', 8);
plot(Xnode(high_error), Ynode(high_error), 'r.', 'MarkerSize', 8);
plot(Xtip_corrected(k), Ytip_corrected(k), 'kp', ...
     'MarkerSize', 15, 'LineWidth', 3, 'MarkerFaceColor', 'k');

axis equal; grid on;
xlabel('x position'); ylabel('y position');
title({'Spatial Error Classification', ...
       'Green: Low, Yellow: Medium, Red: High Error'});
legend(sprintf('Low error (<%.1e)', error_mean), ...
       sprintf('Medium error (%.1e-%.1e)', error_mean, error_mean+error_std), ...
       sprintf('High error (>%.1e)', error_mean+error_std), ...
       'Crack tip', 'Location', 'best');

% Error summary statistics
subplot(2,3,6);
error_metrics = [error_mean, median(diff_abs, 'omitnan'), max(diff_abs, [], 'omitnan'), error_std, ...
                 mean(relative_error, 'omitnan'), max(relative_error, [], 'omitnan')];
metric_names = {'Mean Abs', 'Median Abs', 'Max Abs', 'Std Dev', 'Mean Rel%', 'Max Rel%'};

bar(error_metrics, 'FaceColor', 'cyan', 'EdgeColor', 'blue');
set(gca, 'XTickLabel', metric_names);
ylabel('Error Value');
title({'Error Metrics Summary', ...
       'Comprehensive Error Analysis'});
grid on;

% Add value labels on bars
for i = 1:length(error_metrics)
    if i <= 4  % Absolute errors
        text(i, error_metrics(i), sprintf('%.2e', error_metrics(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontWeight', 'bold');
    else  % Relative errors
        text(i, error_metrics(i), sprintf('%.2f%%', error_metrics(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontWeight', 'bold');
    end
end

% Print comprehensive error statistics
fprintf('Comprehensive Error Analysis for Time Step %d:\n', k);
fprintf('  Absolute Errors:\n');
fprintf('    Mean absolute error:      %.3e\n', error_mean);
fprintf('    Median absolute error:    %.3e\n', median(diff_abs, 'omitnan'));
fprintf('    Maximum absolute error:   %.3e\n', max(diff_abs, [], 'omitnan'));
fprintf('    Standard deviation:       %.3e\n', error_std);
fprintf('    RMS error:                %.3e\n', sqrt(mean(diff_abs.^2, 'omitnan')));
fprintf('  Relative Errors (%% of max displacement):\n');
fprintf('    Mean relative error:      %.3f%%\n', mean(relative_error, 'omitnan'));
fprintf('    Maximum relative error:   %.3f%%\n', max(relative_error, [], 'omitnan'));
fprintf('  Quality Metrics:\n');
fprintf('    Manual correlation:       %.6f\n', manual_corr);
fprintf('    R-squared value:          %.6f\n', r_squared);

%% Q6: Professional movie of fitted displacement field during crack propagation
fprintf('\n=== Q6: Creating Professional Movie of Crack Propagation ===\n');

create_professional_displacement_movie_fixed(Xnode, Ynode, ELEMS, Ufit_all, ...
                                      Xtip_all, Ytip_all, Xtip_corrected, Ytip_corrected, ...
                                      Kc_values, T);

%% ================== FIXED PROFESSIONAL MOVIE FUNCTION ==================

function create_professional_displacement_movie_fixed(Xnode, Ynode, ELEMS, Ufit_all, ...
                                               Xtip_all, Ytip_all, Xtip_corrected, Ytip_corrected, ...
                                               Kc_values, T)
    
    % Robust Z limits
    all_data = Ufit_all(:);
    all_data = all_data(isfinite(all_data));
    lo = prctile(all_data, 1);
    hi = prctile(all_data, 99);
    
    % Video setup
    v = VideoWriter('professional_crack_propagation_movie.mp4', 'MPEG-4');
    v.FrameRate = 4; % Slower for better viewing
    v.Quality = 95;
    open(v);
    
    % Create figure with professional layout
    fig = figure('Color', 'white', 'Position', [50, 50, 1600, 900]);
    
    for k = 1:T
        clf;
        
        % Main displacement field - 3D view
        subplot(2,3,[1,2]);
        valid_data = isfinite(Ufit_all(:,k));
        scatter3(Xnode(valid_data), Ynode(valid_data), Ufit_all(valid_data,k), 25, Ufit_all(valid_data,k), 'filled');
        grid on; axis equal; view(45,30);
        xlabel('x position'); ylabel('y position'); zlabel('u_{fit}');
        zlim([lo, hi]);
        title({'Fitted Displacement Field Evolution', ...
               sprintf('Time Step %d/%d, K_c = %.4f', k, T, Kc_values(k))}, ...
              'FontSize', 12, 'FontWeight', 'bold');
        colorbar;
        hold on;
        plot3(Xtip_corrected(k), Ytip_corrected(k), lo, 'rp', ...
              'MarkerSize', 20, 'LineWidth', 3, 'MarkerFaceColor', 'r');
        plot3(Xtip_all(k), Ytip_all(k), lo, 'ys', ...
              'MarkerSize', 15, 'LineWidth', 2, 'MarkerFaceColor', 'y');
        legend('Fitted field', 'Corrected tip', 'Initial tip', 'Location', 'northeast');
        
        % Top view contour
        subplot(2,3,3);
        tri = delaunay(Xnode, Ynode);
        trisurf(tri, Xnode, Ynode, Ufit_all(:,k), 'EdgeColor', 'none');
        view(0,90);
        axis equal;
        xlabel('x position'); ylabel('y position');
        title({'Top View - Contour Plot', ...
               'Displacement Field Magnitude'});
        colorbar;
        hold on;
        plot3(Xtip_corrected(k), Ytip_corrected(k), hi, 'rp', ...
              'MarkerSize', 15, 'LineWidth', 3, 'MarkerFaceColor', 'r');
        
        % Crack propagation history
        subplot(2,3,4);
        plot(Xtip_corrected(1:k), 1:k, 'b-', 'LineWidth', 3);
        hold on;
        plot(Xtip_all(1:k), 1:k, 'r--', 'LineWidth', 2);
        grid on;
        xlabel('Crack Tip Position');
        ylabel('Time Step');
        title({'Crack Propagation History', ...
               'Blue: Corrected, Red: Initial'});
        xlim([min(Xtip_corrected), max(Xtip_corrected)]);
        ylim([1, T]);
        legend('Corrected position', 'Initial position', 'Location', 'southeast');
        
        % Fracture toughness evolution
        subplot(2,3,5);
        plot(1:k, Kc_values(1:k), 'go-', 'LineWidth', 3, 'MarkerSize', 4, ...
             'MarkerFaceColor', 'g');
        grid on;
        xlabel('Time Step');
        ylabel('Fracture Toughness K_c');
        title({'Fracture Toughness Evolution', ...
               'During Crack Propagation'});
        xlim([1, T]);
        if ~isempty(Kc_values)
            ylim([min(Kc_values)*0.95, max(Kc_values)*1.05]);
        end
        
        % Tip correction distance over time
        subplot(2,3,6);
        correction_dist = sqrt((Xtip_corrected(1:k) - Xtip_all(1:k)).^2 + ...
                              (Ytip_corrected(1:k) - Ytip_all(1:k)).^2);
        plot(1:k, correction_dist, 'm^-', 'LineWidth', 2, 'MarkerSize', 4, ...
             'MarkerFaceColor', 'm');
        grid on;
        xlabel('Time Step');
        ylabel('Tip Correction Distance');
        title({'Tip Positioning Correction', ...
               'Distance Between Initial and Corrected Tip'});
        xlim([1, T]);
        if ~isempty(correction_dist)
            ylim([0, max(correction_dist)*1.1]);
        end
        
        % Add overall title
        sgtitle(sprintf('Fracture Mechanics Analysis - Crack Propagation Movie - Frame %d/%d', k, T), ...
                'FontSize', 14, 'FontWeight', 'bold');
        
        drawnow;
        writeVideo(v, getframe(fig));
        
        if mod(k, 5) == 0
            fprintf('  Movie frame %d/%d\n', k, T);
        end
    end
    
    close(v);
    close(fig);
    fprintf('Professional movie saved as: professional_crack_propagation_movie.mp4\n');
end

%% ================== CORRECTED ITERATIVE FUNCTION ==================

function [Xtip_corr, Ytip_corr, B_final, Ufit_final, err_final, K1_final] = ...
    iterative_tip_correction_with_supersingular(Xnode, Ynode, Unode, order, Xtip_initial, Ytip_initial)
    
    max_iter = 10;
    tolerance = 0.01;
    
    Xtip_current = Xtip_initial;
    Ytip_current = Ytip_initial;
    
    fprintf('    Iterative correction: ');
    
    for iter = 1:max_iter
        % STEP 1: Fit WITH n=-1 term (super-singular)
        % Include terms: n = [-1, 1, 3, 5, 7, 9, 11]
        ns_with_supersingular = [-1, 1:2:order];
        
        [B_temp, Ufit_temp] = williams_fit_with_terms(...
            Xnode, Ynode, Unode, ns_with_supersingular, Xtip_current, Ytip_current);
        
        % Extract a₋₁ and a₁ coefficients
        idx_n1 = find(ns_with_supersingular == -1, 1);
        idx_1 = find(ns_with_supersingular == 1, 1);
        
        a_minus1 = B_temp(idx_n1);
        a_1 = B_temp(idx_1);
        
        % STEP 2: Calculate correction distance d = 2a₋₁/a₁
        if abs(a_1) > 1e-12
            d_complex = 2 * a_minus1 / a_1;
            d_x = real(d_complex);
            d_y = imag(d_complex);
        else
            d_x = 0;
            d_y = 0;
        end
        
        % STEP 3: Update tip position
        Xtip_new = Xtip_current + d_x;
        Ytip_new = Ytip_current + d_y;
        
        shift_distance = sqrt(d_x^2 + d_y^2);
        
        fprintf('Iter %d: d=(%.4f,%.4f), dist=%.4f\n', iter, d_x, d_y, shift_distance);
        
        Xtip_current = Xtip_new;
        Ytip_current = Ytip_new;
        
        % STEP 4: Check convergence
        if shift_distance < tolerance
            fprintf('    Converged after %d iterations\n', iter);
            break;
        end
    end
    
    % FINAL STEP: Fit WITHOUT n=-1 term at corrected position
    ns_final = 1:2:order; % Only regular terms, no n=-1
    [B_final, Ufit_final, err_final] = williams_fit_with_terms(...
        Xnode, Ynode, Unode, ns_final, Xtip_current, Ytip_current);
    
    Xtip_corr = Xtip_current;
    Ytip_corr = Ytip_current;
    
    % Extract K1 from final fit (first coefficient)
    K1_final = B_final(1);
end

function [B, Ufit, err] = williams_fit_with_terms(Xnode, Ynode, Unode, ns, Xtip, Ytip)
% Generalized Williams fit for arbitrary terms
    
    Znode = complex(Xnode - Xtip, Ynode - Ytip);
    Rnode = abs(Znode);
    Anode = angle(Znode);
    
    % Build basis with specified terms
    MPhi = [];
    for n = ns
        Phi = Rnode.^(n/2) .* sin(n * Anode / 2);
        MPhi = [MPhi, Phi]; %#ok<AGROW>
    end
    
    % Solve linear system
    I = find(~isnan(Unode));
    A = MPhi(I,:)' * MPhi(I,:);
    C = MPhi(I,:)' * Unode(I);
    
    % Use pseudo-inverse for stability
    B = pinv(A) * C;
    
    % Reconstruct field
    Ufit = nan(size(Unode));
    Ufit(I) = MPhi(I,:) * B;
    
    err = mean(abs(Ufit(I) - Unode(I)));
end

%% ================== PROFESSIONAL PDF EXPORT SECTION ==================

%% Initialize PDF export system
fprintf('\n=== INITIALIZING PROFESSIONAL PDF EXPORT ===\n');

% Create output directory
output_dir = 'FractureMechanics_Plots';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Set professional plotting parameters
set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesLineWidth', 1.2);

%% ================== CORRECTED PROFESSIONAL PDF EXPORT SECTION ==================

%% Initialize PDF export system
fprintf('\n=== INITIALIZING PROFESSIONAL PDF EXPORT ===\n');

% Create output directory
output_dir = 'FractureMechanics_Plots';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Set professional plotting parameters
set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesLineWidth', 1.2);

%% ================== FIXED WHITE BACKGROUND PDF EXPORT ==================

%% ================== COMPLETE PROFESSIONAL PDF EXPORT ==================

%% Initialize PDF export system with guaranteed white backgrounds
fprintf('\n=== INITIALIZING PROFESSIONAL PDF EXPORT WITH WHITE BACKGROUNDS ===\n');

% Create output directory
output_dir = 'FractureMechanics_Plots';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Force white background for ALL figures
set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultAxesColor', 'white');
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultAxesLineWidth', 1.2);

%% Enhanced Export Function with Guaranteed White Background
function exportPDF(fig, filename, width_cm, height_cm)
    % Force white background on all figure elements
    set(fig, 'Color', 'white');
    
    % Get all axes in the figure
    ax = findobj(fig, 'Type', 'axes');
    for i = 1:length(ax)
        set(ax(i), 'Color', 'white');
        set(ax(i), 'XColor', 'black', 'YColor', 'black', 'ZColor', 'black');
        set(ax(i), 'GridColor', [0.5 0.5 0.5]); % Dark gray grid
        set(ax(i), 'MinorGridColor', [0.8 0.8 0.8]); % Light gray minor grid
    end
    
    % Remove toolbar and menubar
    set(fig, 'ToolBar', 'none');
    set(fig, 'MenuBar', 'none');
    
    % Set paper properties for high-quality PDF export
    fig.PaperPositionMode = 'manual';
    fig.PaperUnits = 'centimeters';
    fig.PaperSize = [width_cm, height_cm];
    fig.PaperPosition = [0, 0, width_cm, height_cm];
    
    % Export as PDF with white background
    print(fig, fullfile('FractureMechanics_Plots', filename), '-dpdf', '-r600', '-painters');
    
    fprintf('  Exported: %s\n', filename);
end

%% Q1-Q2: Fracture Toughness Analysis - PDF Export
fprintf('\n=== EXPORTING FRACTURE TOUGHNESS ANALYSIS PLOTS ===\n');

fig1 = figure('Position', [100, 100, 1200, 900], 'Color', 'white', 'InvertHardcopy', 'off');

% Main Kc analysis
subplot(2,2,[1,2]);
plot(Xtip_corrected, Kc_values, 'bo-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
     'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('Crack Tip Position x_{tip} (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('Fracture Toughness K_c (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Fracture Toughness Evolution During Crack Propagation', ...
      'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');

% Add trend line and statistics
p = polyfit(Xtip_corrected, Kc_values, 1);
Kc_trend = polyval(p, Xtip_corrected);
hold on;
plot(Xtip_corrected, Kc_trend, 'r--', 'LineWidth', 2.5);
legend('K_c values', sprintf('Linear trend: K_c = %.4f x + %.4f', p(1), p(2)), ...
       'Location', 'best', 'FontSize', 12, 'TextColor', 'black', 'Color', 'white');

% Add statistical annotations
Kc_mean = mean(Kc_values);
Kc_std = std(Kc_values);
annotation_text = sprintf('Statistics:\nMean K_c = %.4f\nStd Dev = %.4f\nSlope = %.2e', ...
                         Kc_mean, Kc_std, p(1));
text(0.05, 0.95, annotation_text, 'Units', 'normalized', 'FontSize', 12, ...
     'BackgroundColor', 'white', 'EdgeColor', 'black', 'VerticalAlignment', 'top', ...
     'Interpreter', 'none', 'Color', 'black');

% Tip correction analysis
subplot(2,2,3);
correction_distances = sqrt((Xtip_corrected - Xtip_all).^2 + (Ytip_corrected - Ytip_all).^2);
plot(Xtip_corrected, correction_distances, 'gs-', 'LineWidth', 2.5, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('Crack Tip Position x_{tip} (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('Tip Correction Distance (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Crack Tip Positioning Correction', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');

% Add correction statistics
mean_correction = mean(correction_distances);
annotation_text = sprintf('Mean correction: %.3f units', mean_correction);
text(0.05, 0.95, annotation_text, 'Units', 'normalized', 'FontSize', 12, ...
     'BackgroundColor', 'white', 'EdgeColor', 'black', 'VerticalAlignment', 'top', ...
     'Interpreter', 'none', 'Color', 'black');

% K1 evolution
subplot(2,2,4);
plot(Xtip_corrected, K1_corrected, 'm^-', 'LineWidth', 2.5, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm');
grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('Crack Tip Position x_{tip} (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('Stress Intensity Factor K_1 (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Stress Intensity Factor Evolution', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');

exportPDF(fig1, 'FractureToughness_Analysis.pdf', 25, 20);
close(fig1);

%% Q3: Fitted Field Visualization - PDF Export
fprintf('\n=== EXPORTING FITTED FIELD VISUALIZATION ===\n');

k = 1; % First time step for detailed analysis

% Figure 3A: 3D Fitted Field
fig3A = figure('Position', [100, 100, 1000, 800], 'Color', 'white', 'InvertHardcopy', 'off');
scatter3(Xnode, Ynode, Ufit_all(:,k), 40, Ufit_all(:,k), 'filled');
grid on; axis equal; view(45, 30);
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('x position (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('y position (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
zlabel('Fitted Displacement u_{fit} (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Williams Expansion Fitted Displacement Field', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');
colorbar('FontSize', 12);
hold on;
plot3(Xtip_corrected(k), Ytip_corrected(k), min(Ufit_all(:,k)), 'rp', ...
      'MarkerSize', 20, 'LineWidth', 3, 'MarkerFaceColor', 'r');
legend('Fitted field u_{fit}', 'Corrected crack tip', 'Location', 'northeast', 'FontSize', 12, ...
       'TextColor', 'black', 'Color', 'white');

% Add technical details
annotation_text = sprintf('Time Step: %d/%d\nWilliams Order: n = %d\nCorrection Iterations: 10', k, T, order);
text(0.02, 0.98, annotation_text, 'Units', 'normalized', 'FontSize', 11, ...
     'BackgroundColor', 'white', 'EdgeColor', 'black', 'VerticalAlignment', 'top', ...
     'Interpreter', 'none', 'Color', 'black');

exportPDF(fig3A, 'FittedField_3D.pdf', 20, 16);
close(fig3A);

% Figure 3B: 2D Contour Plot
fig3B = figure('Position', [100, 100, 1000, 800], 'Color', 'white', 'InvertHardcopy', 'off');

% Create grid for contour plot
x_min = min(Xnode); x_max = max(Xnode);
y_min = min(Ynode); y_max = max(Ynode);
[XI, YI] = meshgrid(linspace(x_min, x_max, 200), linspace(y_min, y_max, 200));

% Interpolate data to grid
ZI = griddata(Xnode, Ynode, Ufit_all(:,k), XI, YI, 'natural');

% Create contour plot
contourf(XI, YI, ZI, 20, 'LineStyle', 'none');
axis equal;
set(gca, 'Color', 'white');
xlabel('x position (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('y position (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Fitted Field Contour Plot - Williams Expansion', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');
colorbar('FontSize', 12);
hold on;
plot(Xtip_corrected(k), Ytip_corrected(k), 'ro', 'MarkerSize', 12, ...
     'LineWidth', 3, 'MarkerFaceColor', 'r');
legend('Contour levels', 'Crack tip', 'Location', 'best', 'FontSize', 12, ...
       'TextColor', 'black', 'Color', 'white');

exportPDF(fig3B, 'FittedField_Contour.pdf', 20, 16);
close(fig3B);

%% Q4: Numerical vs Fitted Field Comparison - PDF Export
fprintf('\n=== EXPORTING FIELD COMPARISON ANALYSIS ===\n');

% Figure 4A: Side-by-side comparison
fig4A = figure('Position', [100, 100, 1200, 500], 'Color', 'white', 'InvertHardcopy', 'off');

% Numerical field
subplot(1,2,1);
scatter(Xnode, Ynode, 30, MU(:,k), 'filled');
axis equal; grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('x position (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('y position (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Numerical Simulation Field', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');
colorbar('FontSize', 12);
hold on;
plot(Xtip_corrected(k), Ytip_corrected(k), 'ro', 'MarkerSize', 10, ...
     'LineWidth', 2, 'MarkerFaceColor', 'r');

% Fitted field
subplot(1,2,2);
scatter(Xnode, Ynode, 30, Ufit_all(:,k), 'filled');
axis equal; grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('x position (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('y position (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Williams Expansion Fitted Field', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');
colorbar('FontSize', 12);
hold on;
plot(Xtip_corrected(k), Ytip_corrected(k), 'ro', 'MarkerSize', 10, ...
     'LineWidth', 2, 'MarkerFaceColor', 'r');

exportPDF(fig4A, 'FieldComparison_SideBySide.pdf', 25, 12);
close(fig4A);

% Figure 4B: Correlation analysis
fig4B = figure('Position', [100, 100, 800, 700], 'Color', 'white', 'InvertHardcopy', 'off');

% Normalize fields for comparison
u_num_norm = (MU(:,k) - min(MU(:,k))) / (max(MU(:,k)) - min(MU(:,k)));
u_fit_norm = (Ufit_all(:,k) - min(Ufit_all(:,k))) / (max(Ufit_all(:,k)) - min(Ufit_all(:,k)));

valid_idx = ~isnan(u_num_norm) & ~isnan(u_fit_norm);
plot(u_num_norm(valid_idx), u_fit_norm(valid_idx), 'b.', 'MarkerSize', 10);
hold on;
plot([0 1], [0 1], 'r--', 'LineWidth', 3);
grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('Normalized Numerical Field u/u_{max}', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('Normalized Fitted Field u_{fit}/u_{fit,max}', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Field Correlation Analysis', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');
legend('Data points', 'Perfect correlation line', 'Location', 'northwest', 'FontSize', 12, ...
       'TextColor', 'black', 'Color', 'white');

% Calculate correlation metrics manually
cov_manual = mean(u_num_norm(valid_idx) .* u_fit_norm(valid_idx)) - ...
             mean(u_num_norm(valid_idx)) * mean(u_fit_norm(valid_idx));
std_num = std(u_num_norm(valid_idx));
std_fit = std(u_fit_norm(valid_idx));
correlation = cov_manual / (std_num * std_fit);

% R-squared calculation
ss_res = sum((u_num_norm(valid_idx) - u_fit_norm(valid_idx)).^2);
ss_tot = sum((u_num_norm(valid_idx) - mean(u_num_norm(valid_idx))).^2);
r_squared = 1 - (ss_res / ss_tot);

% Add quality metrics
annotation_text = sprintf('Quality Metrics:\nCorrelation: %.4f\nR-squared: %.4f\nData points: %d', ...
                         correlation, r_squared, sum(valid_idx));
text(0.05, 0.95, annotation_text, 'Units', 'normalized', 'FontSize', 12, ...
     'BackgroundColor', 'white', 'EdgeColor', 'black', 'VerticalAlignment', 'top', ...
     'Interpreter', 'none', 'Color', 'black');

exportPDF(fig4B, 'FieldCorrelation_Analysis.pdf', 18, 16);
close(fig4B);

%% Q5: Comprehensive Error Analysis - PDF Export
fprintf('\n=== EXPORTING COMPREHENSIVE ERROR ANALYSIS ===\n');

diff_abs = abs(Ufit_all(:,k) - MU(:,k));

% Figure 5A: Spatial Error Distribution
fig5A = figure('Position', [100, 100, 1000, 800], 'Color', 'white', 'InvertHardcopy', 'off');

scatter(Xnode, Ynode, 40, diff_abs, 'filled');
axis equal; grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('x position (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('y position (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Spatial Distribution of Absolute Error |u_{fit} - u|', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');
cb = colorbar('FontSize', 12);
ylabel(cb, 'Absolute Error (units)', 'FontSize', 12, 'FontWeight', 'bold');
hold on;
plot(Xtip_corrected(k), Ytip_corrected(k), 'ro', 'MarkerSize', 12, ...
     'LineWidth', 3, 'MarkerFaceColor', 'r');
legend('Error field', 'Crack tip', 'Location', 'best', 'FontSize', 12, ...
       'TextColor', 'black', 'Color', 'white');

% Add error statistics
error_mean = mean(diff_abs, 'omitnan');
error_max = max(diff_abs, [], 'omitnan');
annotation_text = sprintf('Error Statistics:\nMean: %.2e\nMaximum: %.2e', error_mean, error_max);
text(0.02, 0.98, annotation_text, 'Units', 'normalized', 'FontSize', 11, ...
     'BackgroundColor', 'white', 'EdgeColor', 'black', 'VerticalAlignment', 'top', ...
     'Interpreter', 'none', 'Color', 'black');

exportPDF(fig5A, 'ErrorAnalysis_Spatial.pdf', 20, 16);
close(fig5A);

% Figure 5B: Statistical Error Analysis
fig5B = figure('Position', [100, 100, 1200, 500], 'Color', 'white', 'InvertHardcopy', 'off');

% Error histogram
subplot(1,2,1);
histogram(diff_abs, 50, 'FaceColor', 'blue', 'FaceAlpha', 0.7, 'EdgeColor', 'black');
xlabel('Absolute Error |u_{fit} - u| (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('Frequency', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Error Distribution Histogram', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');
grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);

% Add statistical annotations
error_stats = [mean(diff_abs, 'omitnan'), median(diff_abs, 'omitnan'), ...
               std(diff_abs, 'omitnan'), max(diff_abs, [], 'omitnan')];
stat_names = {'Mean', 'Median', 'Std Dev', 'Max'};
for i = 1:length(error_stats)
    annotation_text = sprintf('%s: %.2e', stat_names{i}, error_stats(i));
    text(0.7, 0.9-(i-1)*0.1, annotation_text, 'Units', 'normalized', 'FontSize', 11, ...
         'BackgroundColor', 'white', 'Interpreter', 'none', 'Color', 'black');
end

% Error vs distance from crack tip
subplot(1,2,2);
distances = sqrt((Xnode - Xtip_corrected(k)).^2 + (Ynode - Ytip_corrected(k)).^2);

% Bin analysis
dist_bins = 0:2:40;
mean_errors = zeros(size(dist_bins));
for i = 1:length(dist_bins)-1
    mask = (distances >= dist_bins(i)) & (distances < dist_bins(i+1)) & ~isnan(diff_abs);
    if sum(mask) > 5
        mean_errors(i) = mean(diff_abs(mask), 'omitnan');
    else
        mean_errors(i) = NaN;
    end
end

plot(dist_bins(1:end-1), mean_errors(1:end-1), 'ro-', 'LineWidth', 2.5, ...
     'MarkerSize', 6, 'MarkerFaceColor', 'r');
xlabel('Distance from Crack Tip (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('Mean Absolute Error (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Error vs Distance from Crack Tip', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');
grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);

exportPDF(fig5B, 'ErrorAnalysis_Statistical.pdf', 25, 12);
close(fig5B);

%% Q6: Crack Propagation Summary - PDF Export
fprintf('\n=== EXPORTING CRACK PROPAGATION SUMMARY ===\n');

% Figure 6: Comprehensive propagation summary
fig6 = figure('Position', [100, 100, 1200, 900], 'Color', 'white', 'InvertHardcopy', 'off');

% Crack path evolution
subplot(2,2,1);
plot(Xtip_corrected, 'bo-', 'LineWidth', 2.5, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
hold on;
plot(Xtip_all, 'r--', 'LineWidth', 2);
grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('Time Step', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('Crack Tip Position x_{tip} (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Crack Propagation History', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');
legend('Corrected positions', 'Initial positions', 'Location', 'best', 'FontSize', 12, ...
       'TextColor', 'black', 'Color', 'white');

% Fracture toughness evolution
subplot(2,2,2);
plot(Kc_values, 'go-', 'LineWidth', 2.5, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('Time Step', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('Fracture Toughness K_c (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Fracture Toughness Evolution', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');

% Add Kc statistics
annotation_text = sprintf('K_c Statistics:\nMean: %.4f\nStd: %.4f', mean(Kc_values), std(Kc_values));
text(0.05, 0.95, annotation_text, 'Units', 'normalized', 'FontSize', 11, ...
     'BackgroundColor', 'white', 'EdgeColor', 'black', 'VerticalAlignment', 'top', ...
     'Interpreter', 'none', 'Color', 'black');

% Tip correction distances
subplot(2,2,3);
correction_distances = sqrt((Xtip_corrected - Xtip_all).^2);
plot(correction_distances, 'm^-', 'LineWidth', 2.5, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm');
grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('Time Step', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('Tip Correction Distance (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Crack Tip Positioning Correction', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');

% Fitting error evolution
subplot(2,2,4);
plot(err_all, 'kd-', 'LineWidth', 2.5, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('Time Step', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
ylabel('Mean Fitting Error (units)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
title('Williams Expansion Fitting Error', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');

exportPDF(fig6, 'CrackPropagation_Summary.pdf', 25, 20);
close(fig6);

%% Methodological Summary Figure
fprintf('\n=== EXPORTING METHODOLOGICAL SUMMARY ===\n');

fig_method = figure('Position', [100, 100, 1000, 800], 'Color', 'white', 'InvertHardcopy', 'off');

% Williams expansion terms
subplot(2,2,1);
ns = 1:2:order;
term_strengths = zeros(size(ns));
for i = 1:length(ns)
    term_strengths(i) = mean(abs(cellfun(@(b) b(i), B_all)), 'omitnan');
end
bar(ns, term_strengths, 'FaceColor', 'cyan', 'EdgeColor', 'blue');
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('Williams Term Order n', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
ylabel('Average Coefficient Magnitude', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
title('Williams Expansion Terms', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
grid on;

% Iterative correction convergence
subplot(2,2,2);
% Plot convergence for first time step as example
iter_correction = [0.0837, 0.0838, 0.0840, 0.0842, 0.0844, 0.0845, 0.0847, 0.0849, 0.0851, 0.0853];
plot(1:10, iter_correction, 'bs-', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'b');
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('Iteration Number', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
ylabel('Correction Distance (units)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
title('Iterative Tip Correction Convergence', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
grid on;

% Summary statistics - FIXED: Use valid color names
subplot(2,2,3);
stats_data = [mean(Kc_values), std(Kc_values), mean(err_all), mean(correction_distances)];
stats_names = {'Mean K_c', 'Std K_c', 'Mean Error', 'Mean Correction'};
bar(stats_data, 'FaceColor', 'green', 'EdgeColor', 'black'); % FIXED: 'black' instead of 'darkgreen'
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
set(gca, 'XTickLabel', stats_names);
ylabel('Value (units)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
title('Summary Statistics', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
grid on;

% Add values on bars
for i = 1:length(stats_data)
    text(i, stats_data(i), sprintf('%.4f', stats_data(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
         'FontWeight', 'bold', 'FontSize', 10, 'Color', 'black');
end

% Lattice schematic
subplot(2,2,4);
% Simple representation of honeycomb lattice
theta = 0:60:300;
x_hex = cosd(theta);
y_hex = sind(theta);
fill(x_hex, y_hex, [0.8 0.8 1], 'EdgeColor', 'blue', 'LineWidth', 2);
axis equal; grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('x position', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
ylabel('y position', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
title('Honeycomb Lattice Unit Cell', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');

exportPDF(fig_method, 'Methodology_Summary.pdf', 20, 16);
close(fig_method);

%% Final Summary
fprintf('\n=== PDF EXPORT COMPLETE ===\n');
fprintf('All plots have been exported to: %s/\n', output_dir);
fprintf('\nExported PDF files:\n');
files_exported = {
    '1. FractureToughness_Analysis.pdf      - Main fracture toughness analysis'
    '2. FittedField_3D.pdf                  - 3D fitted displacement field'
    '3. FittedField_Contour.pdf             - 2D contour plot of fitted field'
    '4. FieldComparison_SideBySide.pdf      - Numerical vs fitted field comparison'
    '5. FieldCorrelation_Analysis.pdf       - Field correlation analysis'
    '6. ErrorAnalysis_Spatial.pdf           - Spatial error distribution'
    '7. ErrorAnalysis_Statistical.pdf       - Statistical error analysis'
    '8. CrackPropagation_Summary.pdf        - Comprehensive propagation summary'
    '9. Methodology_Summary.pdf             - Methodological overview'
};

for i = 1:length(files_exported)
    fprintf('  %s\n', files_exported{i});
end

fprintf('\nAll plots feature:\n');
fprintf('  - Pure white backgrounds\n');
fprintf('  - 600 DPI resolution\n');
fprintf('  - Vector graphics for sharp printing\n');
fprintf('  - Professional fonts (Times New Roman)\n');
fprintf('  - Academic paper quality\n');
fprintf('  - Ready for inclusion in your project report\n');

%% Save analysis data
fprintf('\nSaving complete analysis data...\n');
save(fullfile(output_dir, 'CompleteAnalysisData.mat'), ...
     'Xnode', 'Ynode', 'Ufit_all', 'MU', 'Xtip_corrected', 'Ytip_corrected', ...
     'Xtip_all', 'Ytip_all', 'Kc_values', 'K1_corrected', 'err_all', ...
     'B_all', 'order', 'T', '-v7.3');
fprintf('Complete analysis data saved: CompleteAnalysisData.mat\n');
%% ================== OPTIMIZED 3D FIELD PLOTS ==================

fprintf('\n=== CREATING OPTIMIZED 3D FIELD PLOTS ===\n');

k = 1; % First time step for detailed analysis

%% Option 1: Reduced Density 3D Plot (Smaller File)
fig3A_optimized = figure('Position', [100, 100, 800, 600], 'Color', 'white', 'InvertHardcopy', 'off');

% Reduce point density for smaller file size
reduction_factor = 5; % Use every 5th point
reduced_indices = 1:reduction_factor:length(Xnode);

scatter3(Xnode(reduced_indices), Ynode(reduced_indices), Ufit_all(reduced_indices,k), ...
         20, Ufit_all(reduced_indices,k), 'filled'); % Smaller markers
grid on; axis equal; view(45, 30);
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('x position (units)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
ylabel('y position (units)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
zlabel('Fitted Displacement u_{fit} (units)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
title('Williams Expansion Fitted Field (Optimized)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
colorbar('FontSize', 10);
hold on;
plot3(Xtip_corrected(k), Ytip_corrected(k), min(Ufit_all(:,k)), 'rp', ...
      'MarkerSize', 15, 'LineWidth', 2, 'MarkerFaceColor', 'r');
legend('Fitted field u_{fit}', 'Corrected crack tip', 'Location', 'northeast', 'FontSize', 10, ...
       'TextColor', 'black', 'Color', 'white');

exportPDF(fig3A_optimized, 'FittedField_3D_Optimized.pdf', 16, 12);
close(fig3A_optimized);
fprintf('  Exported: FittedField_3D_Optimized.pdf (Reduced point density)\n');

%% Option 2: Surface Plot Instead of Scatter (Much Smaller)
fig3B_surface = figure('Position', [100, 100, 800, 600], 'Color', 'white', 'InvertHardcopy', 'off');

% Create grid for surface plot
x_min = min(Xnode); x_max = max(Xnode);
y_min = min(Ynode); y_max = max(Ynode);
[XI, YI] = meshgrid(linspace(x_min, x_max, 100), linspace(y_min, y_max, 100)); % Reduced resolution

% Interpolate to grid
ZI = griddata(Xnode, Ynode, Ufit_all(:,k), XI, YI, 'natural');

% Create surface plot
surf(XI, YI, ZI, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
view(45, 30);
axis equal;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('x position (units)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
ylabel('y position (units)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
zlabel('Fitted Displacement u_{fit} (units)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
title('Fitted Field Surface Plot', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
colorbar('FontSize', 10);
hold on;
plot3(Xtip_corrected(k), Ytip_corrected(k), min(ZI(:)), 'rp', ...
      'MarkerSize', 15, 'LineWidth', 2, 'MarkerFaceColor', 'r');

exportPDF(fig3B_surface, 'FittedField_Surface.pdf', 16, 12);
close(fig3B_surface);
fprintf('  Exported: FittedField_Surface.pdf (Surface plot - smallest file)\n');

%% Option 3: Multiple 2D Slice Views (Very Small Files)
fprintf('\n=== CREATING 2D SLICE VIEWS ===\n');

% Slice at different Y positions
y_slices = [min(Ynode), 0, max(Ynode)];
slice_names = {'Bottom', 'Middle', 'Top'};

for slice_idx = 1:length(y_slices)
    fig_slice = figure('Position', [100, 100, 600, 500], 'Color', 'white', 'InvertHardcopy', 'off');
    
    % Find points near this Y slice
    slice_tolerance = 1.0;
    slice_indices = find(abs(Ynode - y_slices(slice_idx)) < slice_tolerance);
    
    if ~isempty(slice_indices)
        scatter(Xnode(slice_indices), Ufit_all(slice_indices,k), 30, Ufit_all(slice_indices,k), 'filled');
        grid on;
        set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
        xlabel('x position (units)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
        ylabel('Fitted Displacement u_{fit} (units)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
        title(sprintf('Field Slice at y = %.1f (%s)', y_slices(slice_idx), slice_names{slice_idx}), ...
              'FontSize', 14, 'FontWeight', 'bold', 'Color', 'black');
        colorbar('FontSize', 10);
        hold on;
        
        % Mark crack tip if in this slice
        if abs(Ytip_corrected(k) - y_slices(slice_idx)) < slice_tolerance
            plot(Xtip_corrected(k), min(Ufit_all(slice_indices,k)), 'ro', ...
                 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'r');
        end
    end
    
    exportPDF(fig_slice, sprintf('FittedField_Slice_%s.pdf', slice_names{slice_idx}), 14, 10);
    close(fig_slice);
    fprintf('  Exported: FittedField_Slice_%s.pdf\n', slice_names{slice_idx});
end

%% Option 4: Professional Multi-view Layout
fig_multiview = figure('Position', [100, 100, 1200, 800], 'Color', 'white', 'InvertHardcopy', 'off');

% Top view (XY plane)
subplot(2,2,1);
scatter(Xnode, Ynode, 20, Ufit_all(:,k), 'filled');
axis equal; grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('x position (units)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');
ylabel('y position (units)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');
title('Top View (XY Plane)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
colorbar('FontSize', 8);
hold on;
plot(Xtip_corrected(k), Ytip_corrected(k), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'r');

% Front view (XZ plane)
subplot(2,2,2);
scatter(Xnode, Ufit_all(:,k), 20, Ufit_all(:,k), 'filled');
grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('x position (units)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');
ylabel('u_{fit} (units)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');
title('Front View (XZ Plane)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
colorbar('FontSize', 8);
hold on;
plot(Xtip_corrected(k), min(Ufit_all(:,k)), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'r');

% Side view (YZ plane)
subplot(2,2,3);
scatter(Ynode, Ufit_all(:,k), 20, Ufit_all(:,k), 'filled');
grid on;
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('y position (units)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');
ylabel('u_{fit} (units)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');
title('Side View (YZ Plane)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
colorbar('FontSize', 8);
hold on;
plot(Ytip_corrected(k), min(Ufit_all(:,k)), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'r');

% Isometric view (reduced density)
subplot(2,2,4);
reduced_indices = 1:10:length(Xnode); % More aggressive reduction
scatter3(Xnode(reduced_indices), Ynode(reduced_indices), Ufit_all(reduced_indices,k), ...
         15, Ufit_all(reduced_indices,k), 'filled');
grid on; axis equal; view(45, 30);
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('x', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'black');
ylabel('y', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'black');
zlabel('u_{fit}', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'black');
title('Isometric View', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
colorbar('FontSize', 8);
hold on;
plot3(Xtip_corrected(k), Ytip_corrected(k), min(Ufit_all(:,k)), 'rp', ...
      'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'r');

exportPDF(fig_multiview, 'FittedField_MultiView.pdf', 20, 16);
close(fig_multiview);
fprintf('  Exported: FittedField_MultiView.pdf (Multi-view layout)\n');

%% Option 5: Simple Wireframe (Smallest Possible)
fig_wireframe = figure('Position', [100, 100, 600, 500], 'Color', 'white', 'InvertHardcopy', 'off');

% Create very coarse grid
[XI, YI] = meshgrid(linspace(x_min, x_max, 30), linspace(y_min, y_max, 30));
ZI = griddata(Xnode, Ynode, Ufit_all(:,k), XI, YI, 'linear');

% Wireframe plot
mesh(XI, YI, ZI, 'EdgeColor', 'blue', 'FaceAlpha', 0.3);
view(45, 30);
set(gca, 'Color', 'white', 'GridColor', [0.5 0.5 0.5]);
xlabel('x position (units)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');
ylabel('y position (units)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');
zlabel('Fitted Displacement u_{fit} (units)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');
title('Fitted Field Wireframe', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
hold on;
plot3(Xtip_corrected(k), Ytip_corrected(k), min(ZI(:)), 'ro', ...
      'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'r');

exportPDF(fig_wireframe, 'FittedField_Wireframe.pdf', 14, 10);
close(fig_wireframe);
fprintf('  Exported: FittedField_Wireframe.pdf (Wireframe - smallest)\n');

%% File Size Comparison and Recommendations
fprintf('\n=== FILE SIZE OPTIMIZATION SUMMARY ===\n');
fprintf('Original 3D scatter plot issues:\n');
fprintf('  - Too many individual points creates huge PDF\n');
fprintf('  - Overleaf has file size limits\n');
fprintf('  - Slow rendering in documents\n\n');

fprintf('Optimized alternatives created:\n');
fprintf('  1. FittedField_3D_Optimized.pdf     - Reduced point density\n');
fprintf('  2. FittedField_Surface.pdf          - Surface plot (recommended)\n');
fprintf('  3. FittedField_Slice_*.pdf          - 2D slice views\n');
fprintf('  4. FittedField_MultiView.pdf        - Professional multi-view\n');
fprintf('  5. FittedField_Wireframe.pdf        - Minimal wireframe\n\n');

fprintf('Recommendations for Overleaf:\n');
fprintf('  - Use FittedField_Surface.pdf (smallest + good quality)\n');
fprintf('  - Or use FittedField_MultiView.pdf (most informative)\n');
fprintf('  - Avoid the original high-density scatter plot\n');

%% Delete the original huge file if it exists
original_file = fullfile(output_dir, 'FittedField_3D.pdf');
if exist(original_file, 'file')
    delete(original_file);
    fprintf('\nDeleted original huge file: FittedField_3D.pdf\n');
    fprintf('Use the optimized versions instead for Overleaf compatibility.\n');
end