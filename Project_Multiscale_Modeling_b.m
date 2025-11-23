% Project_Multiscale_AMR.m
% ========================================================================
% Fixed Iteration AMR: Uniform vs. Adaptive FEM Comparison (corrected)
% ========================================================================
clear; close all; clc;

%% =====================================================================
% 1. GLOBAL PARAMETERS
% =====================================================================
L       = 1.0;              % Domain size
R0      = 0.01;             % Proton radius (Regularization)
center  = [0.5*L, 0.5*L];   % Proton position
epsilon = 3 * R0;           % Gaussian width
sigma = epsilon / 3;        % std. deviation of Gaussian
norm_factor = 1 / (2*pi*sigma^2); % Normalization factor
rho_func = @(x,y) norm_factor * exp(-((x-center(1)).^2 + (y-center(2)).^2)/(sigma^2));

% Analytical Solution Function (regularized 1/r)
phi_ana_func = @(r) -(1/(2*pi)) * ...
    ( log(r + 1e-30) - 0.5 * expint( r.^2 ./ (2*sigma^2) ) );


%% =====================================================================
% 2. BASELINE UNIFORM MESH SOLUTION (For comparison)
% =====================================================================
fprintf('--- 2. Computing Baseline Uniform Mesh Solution ---\n');
Nx_uniform = 100; % Fine uniform mesh for baseline (reasonable)
Ny_uniform = 100;
[nodes_uni, elem_uni] = generate_initial_mesh(L, Nx_uniform, Ny_uniform);
N_uni = size(nodes_uni, 1);
Ne_uni = size(elem_uni, 1);

[A_uni, b_uni, ~, ~, ~] = assemble_system(nodes_uni, elem_uni, N_uni, Ne_uni, center, rho_func);

% Apply BCs (Dirichlet phi=0)
tol = 1e-12;
boundary_nodes_uni = find( abs(nodes_uni(:,1)) < tol | ...
                           abs(nodes_uni(:,1)-L) < tol | ...
                           abs(nodes_uni(:,2)) < tol | ...
                           abs(nodes_uni(:,2)-L) < tol );

free_nodes_uni = setdiff(1:N_uni, boundary_nodes_uni);

phi_uni = zeros(N_uni, 1);
A_ff_uni = A_uni(free_nodes_uni, free_nodes_uni);
b_f_uni  = b_uni(free_nodes_uni);
phi_uni(free_nodes_uni) = A_ff_uni \ b_f_uni;
phi_uni(boundary_nodes_uni) = 0;  % Dirichlet BC

r_uni = sqrt((nodes_uni(:,1)-center(1)).^2 + (nodes_uni(:,2)-center(2)).^2);
phi_ana_uni = phi_ana_func(r_uni);

%% =====================================================================
% 3. ADAPTIVE MESH REFINEMENT (AMR) SOLUTION (WITH PROGRESS)
% =====================================================================
fprintf('\n--- 3. Computing Adaptive Mesh Refinement (AMR) Solution ---\n');

% AMR Parameters
Nx_init = 50;
Ny_init = 50;
MAX_ITER = 50;
THETA = 0.3;
use_waitbar = true;   % set false if you run in a non-GUI MATLAB (e.g., headless)

[nodes_amr, elem_amr] = generate_initial_mesh(L, Nx_init, Ny_init);

% prepare persistent edge midpoint map (so midpoints are reused)
edge_mid_map = containers.Map('KeyType','char','ValueType','double');

% Initialize fallback containers (store last stable)
prev_nodes = nodes_amr;
prev_elem = elem_amr;
N_amr = size(nodes_amr,1);
phi_amr = zeros(N_amr,1);
prev_phi_amr = phi_amr;
final_iter = 0;

% Initialize waitbar safely
if use_waitbar
    try
        hwait = waitbar(0,'Starting AMR...');
    catch
        hwait = [];
        use_waitbar = false;
    end
else
    hwait = [];
end

for iter = 1:MAX_ITER
    % update waitbar
    if use_waitbar && ~isempty(hwait) && ishandle(hwait)
        waitbar((iter-1)/MAX_ITER, hwait, sprintf('AMR iteration %d/%d', iter, MAX_ITER));
    end

    fprintf('\n[AMR] Iteration %d/%d\n', iter, MAX_ITER);
    N = size(nodes_amr, 1);
    Ne = size(elem_amr, 1);

    %% --------------------- ASSEMBLY ---------------------
    fprintf('[AMR]  Assembling system... ');
    [A, b, dN_dx_list, dN_dy_list, Area_list] = ...
        assemble_system(nodes_amr, elem_amr, N, Ne, center, rho_func);
    fprintf('done.\n');

    %% --------------------- SOLVE ------------------------
    fprintf('[AMR]  Solving linear system... ');
    % Identify Dirichlet boundary nodes correctly (domain [0,L] x [0,L])
    on_boundary = (abs(nodes_amr(:,1)) < tol) | (abs(nodes_amr(:,1)-L) < tol) | ...
                  (abs(nodes_amr(:,2)) < tol) | (abs(nodes_amr(:,2)-L) < tol);
    free_nodes = find(~on_boundary);

    if isempty(free_nodes)
        warning('[AMR] No free nodes available (all nodes fixed). Stopping.');
        break;
    end

    A_ff = A(free_nodes, free_nodes);
    b_f  = b(free_nodes);

    phi_amr = zeros(N,1);          % full-length vector with BC entries = 0
    phi_amr(free_nodes) = A_ff \ b_f;
    fprintf('done.\n');

    %% ------------------ STABILITY CHECK -----------------
    if any(~isfinite(phi_amr))
        fprintf('[AMR]  ⚠ Numerical blow-up detected (NaN/Inf). Stopping and using last stable mesh.\n');
        break;
    end

    % Store stable iteration (fallback)
    prev_nodes = nodes_amr;
    prev_elem  = elem_amr;
    prev_phi_amr = phi_amr;
    final_iter = iter;

    %% ---------------- ERROR ESTIMATION ------------------
    fprintf('[AMR]  Estimating error... ');
    eta = error_estimation(nodes_amr, elem_amr, phi_amr, ...
                           dN_dx_list, dN_dy_list, Area_list, N, Ne);
    Total_Error_L2 = norm(eta);
    fprintf('done. Total Error = %.4e\n', Total_Error_L2);

    %% -------------------- MARK --------------------------
    fprintf('[AMR]  Marking elements... ');
    marked_elements = mark_elements(eta, THETA);
    fprintf('%d elements marked.\n', numel(marked_elements));

    if isempty(marked_elements)
        fprintf('[AMR]  No marked elements -> adaptation complete.\n');
        break;
    end

    %% -------------------- ENFORCE CLOSURE ----------------
    % This ensures conformity: expand marked_elements to their neighbors until closure
    marked_closed = closure_neighbours(marked_elements, elem_amr);
    fprintf('[AMR]  After closure: %d elements to refine.\n', numel(marked_closed));

    %% -------------------- REFINE ------------------------
    fprintf('[AMR]  Refining mesh (old Ne=%d)... ', Ne);
    % NOTE: refine_mesh must accept and return edge_mid_map
    [nodes_amr, elem_amr, edge_mid_map] = refine_mesh(nodes_amr, elem_amr, marked_closed, edge_mid_map);
    fprintf('new Ne=%d.\n', size(elem_amr,1));

    % quick safety
    if size(nodes_amr,1) > 5e5
        warning('[AMR] Node count exceeded safety limit. Stopping.');
        break;
    end
end

% close waitbar if present
if use_waitbar && ~isempty(hwait) && ishandle(hwait)
    close(hwait);
end

% If we broke early due to NaN/Inf, use the last stable result (prev_*)
nodes_final = prev_nodes;
elem_final  = prev_elem;
phi_final   = prev_phi_amr;   % full nodal vector matching nodes_final
N_final = size(nodes_final,1);

fprintf('\nUsing stable solution from AMR iteration %d\n', final_iter);
fprintf('Final AMR mesh: N_nodes = %d, N_elements = %d\n', N_final, size(elem_final,1));

% Safety checks
if max(elem_final(:)) > N_final
    error('FINAL MESH ERROR: element references node %d but only %d nodes exist', max(elem_final(:)), N_final);
end
if length(phi_final) ~= N_final
    error('FINAL MESH MISMATCH: phi_final length (%d) != number of nodes (%d).', length(phi_final), N_final);
end

% Build triangulations
TR_final = triangulation(elem_final, nodes_final);
TR_uni   = triangulation(elem_uni, nodes_uni);

% Define probe line (before interpolation)
line_x = linspace(0, L, 500);
line_y = center(2) * ones(size(line_x));

% Interpolate along probe line (safe)
disp('[AMR] Interpolating phi along probe line (final stable mesh)...');
phi_amr_line = triinterp(TR_final, phi_final, [line_x', line_y']);
phi_uni_line = triinterp(TR_uni, phi_uni, [line_x', line_y']);

% Compute analytical along line
r_line = abs(line_x - center(1));
phi_ana_line = phi_ana_func(r_line);

% Prepare phi_plot for final mesh (zeros on Dirichlet boundary)
tol = 1e-12;
boundary_nodes_final = find( abs(nodes_final(:,1)) < tol | abs(nodes_final(:,1)-L) < tol | ...
                             abs(nodes_final(:,2)) < tol | abs(nodes_final(:,2)-L) < tol );
free_nodes_final = setdiff(1:N_final, boundary_nodes_final);

phi_plot_amr = zeros(N_final,1);
phi_plot_amr(free_nodes_final) = phi_final(free_nodes_final);
phi_plot_amr(boundary_nodes_final) = 0;
phi_plot_amr = min(phi_plot_amr, 5); % clip for visualization (phi_clip_value = 5)

%% =====================================================================
% 4. PLOTTING, ERROR COMPUTATION & ADAPTIVE DIAGNOSTICS
% =====================================================================

% --- Physical constants & scaling (change phys_scale if your domain is different)
eps0       = 8.8541878128e-12;    % vacuum permittivity (F/m)
kCoulomb   = 1/(4*pi*eps0);       % Coulomb constant
phys_scale = 1e-13;               % length in meters corresponding to unit domain (set to 1e-13 m if L=1 -> domain [0,1] maps to [0,1e-13] m)
R0_phys    = 1e-15;               % physical proton "regularization" radius (m)
% note: in your script 'R0' was a normalized regularization (e.g. 0.01). We keep physical R0 separate.

% physical center
center_phys = center * phys_scale;

% Analytical potential (regularized Coulomb) in SI units
phi_analytic_phys = @(x_phys, y_phys) (kCoulomb) ./ sqrt((x_phys - center_phys(1)).^2 + (y_phys - center_phys(2)).^2 + R0_phys^2);

% --- Ensure nodal vectors used for plotting match node arrays
% Use phi_plot_amr (already prepared) for AMR visual height
% Use phi_uni for uniform solution (already computed earlier)

% clamp for display to avoid giant spikes
phi_clip_value_plot = 2;   % display clip (V) — adjust to taste

% --- Compute analytical values on the same node sets (physical mapping) ---
nodes_uni_phys  = nodes_uni  * phys_scale;
nodes_final_phys= nodes_final* phys_scale;   % final stable AMR mesh nodes

phi_exact_uni   = phi_analytic_phys(nodes_uni_phys(:,1), nodes_uni_phys(:,2));
phi_exact_final = phi_analytic_phys(nodes_final_phys(:,1), nodes_final_phys(:,2));

% --- Prepare safe line data (replace NaNs/interpolate if needed) ---
phi_uni_line_plot = phi_uni_line;
phi_uni_line_plot(isnan(phi_uni_line_plot)) = phi_clip_value_plot;
phi_amr_line_plot = phi_amr_line;
phi_amr_line_plot(isnan(phi_amr_line_plot)) = phi_clip_value_plot;

% -----------------------
% 2D / 3D Visualizations
% -----------------------

% AMR final: 2D (top view)
figure('Name','AMR Final (2D)');
trisurf(elem_final, nodes_final(:,1), nodes_final(:,2), min(phi_plot_amr, phi_clip_value_plot), 'EdgeColor','none');
view(2); axis equal tight; colorbar;
title(sprintf('AMR Final Solution (iteration %d) — clipped to %.2g V', final_iter, phi_clip_value_plot));
xlabel('x (normalized)'); ylabel('y (normalized)');

% AMR final: 3D
figure('Name','AMR Final (3D)');
trisurf(elem_final, nodes_final(:,1), nodes_final(:,2), min(phi_plot_amr, phi_clip_value_plot), ...
        'EdgeColor','none', 'FaceColor','interp');
shading interp; colormap(parula); colorbar;
view(35,40); axis tight equal;
title('AMR Final Solution (3D view)');

% Uniform numerical: 3D
figure('Name','Uniform Numerical (3D)');
trisurf(elem_uni, nodes_uni(:,1), nodes_uni(:,2), min(phi_uni, phi_clip_value_plot), ...
        'EdgeColor','none', 'FaceColor','interp');
shading interp; colormap(jet); colorbar;
view(35,40); axis tight equal;
title('Uniform Mesh Numerical Potential (3D)');

% Analytical on the uniform mesh nodes: 3D (for direct visual comparison)
figure('Name','Analytical (3D) - Uniform nodes');
trisurf(elem_uni, nodes_uni(:,1), nodes_uni(:,2), min(phi_exact_uni, phi_clip_value_plot), ...
        'EdgeColor','none', 'FaceColor','interp');
shading interp; colormap(jet); colorbar;
view(35,40); axis tight equal;
title('Analytical Regularized Coulomb (sampled on uniform nodes)');

% Side-by-side numerical comparison (Uniform vs Adaptive)
figure('Name','Uniform vs Adaptive (3D)');
subplot(1,2,1);
trisurf(elem_uni, nodes_uni(:,1), nodes_uni(:,2), min(phi_uni, phi_clip_value_plot), 'EdgeColor','none'); shading interp;
title('Uniform (numerical)'); view(35,40); axis tight equal; colorbar;

subplot(1,2,2);
trisurf(elem_final, nodes_final(:,1), nodes_final(:,2), min(phi_plot_amr, phi_clip_value_plot), 'EdgeColor','none'); shading interp;
title('Adaptive (numerical)'); view(35,40); axis tight equal; colorbar;

% -----------------------
% Adaptive diagnostic: color by nodal abs-error (log10)
% -----------------------
% nodal absolute error on AMR nodes (use phi_plot_amr which has BC zeros)
nodal_abs_err_final = abs(phi_plot_amr - phi_exact_final);
% avoid zero in log
cdata_err = log10(nodal_abs_err_final + eps);

figure('Name','AMR solution (height) with color = log10(|error|)');
trisurf(elem_final, nodes_final(:,1), nodes_final(:,2), min(phi_plot_amr, phi_clip_value_plot), cdata_err, ...
        'EdgeColor','none', 'FaceColor','interp');
shading interp; colorbar;
colormap(jet);
title('AMR height = \phi_{num}, color = log10(|\phi_{num}-\phi_{ana}|)');
view(35,40); axis tight equal;

%% ---------------- Cross-section with nodes and error ----------------

% Define a cross-section line (horizontal through the proton center)
line_x = linspace(0, L, 500);
line_y = center(2) * ones(size(line_x));

% Interpolate phi along the line
phi_uni_line = triinterp(TR_uni, phi_uni, [line_x', line_y']);
phi_amr_line = triinterp(TR_final, phi_final, [line_x', line_y']);
phi_ana_line = phi_ana_func(abs(line_x - center(1)));

% Replace NaNs for plotting
phi_uni_line_plot = phi_uni_line; phi_uni_line_plot(isnan(phi_uni_line_plot)) = 0;
phi_amr_line_plot = phi_amr_line; phi_amr_line_plot(isnan(phi_amr_line_plot)) = 0;

% Find nodes close to the cross-section line for markers
tol_line = 1e-5;
nodes_on_line_uni = abs(nodes_uni(:,2) - center(2)) < tol_line;
nodes_on_line_amr = abs(nodes_final(:,2) - center(2)) < tol_line;

% Compute error along the line
err_uni_line = abs(phi_uni_line_plot - phi_ana_line);
err_amr_line = abs(phi_amr_line_plot - phi_ana_line);

% ---------------- Plot potential and nodes ----------------
figure('Name','Cross-section with nodes and error');

subplot(2,1,1);
plot(line_x, phi_ana_line, 'r:', 'LineWidth', 2); hold on;
plot(line_x, phi_uni_line_plot, 'b-', 'LineWidth', 2);
plot(line_x, phi_amr_line_plot, 'g--', 'LineWidth', 2);
% show nodes as markers
plot(nodes_uni(nodes_on_line_uni,1), phi_uni(nodes_on_line_uni), 'bo', 'MarkerSize',4, 'MarkerFaceColor','b');
plot(nodes_final(nodes_on_line_amr,1), phi_final(nodes_on_line_amr), 'go', 'MarkerSize',4, 'MarkerFaceColor','g');
xlabel('x [m]'); ylabel('\phi'); title('Cross-section Potential'); grid on;
legend('Analytical', 'Uniform Mesh', 'Adaptive Mesh', 'Uniform Nodes', 'AMR Nodes', 'Location','best');

% ---------------- Plot error ----------------
subplot(2,1,2);
plot(line_x, err_uni_line, 'b-', 'LineWidth', 2); hold on;
plot(line_x, err_amr_line, 'g--', 'LineWidth', 2);
xlabel('x [m]'); ylabel('Error |phi - phi_{ana}|'); title('Cross-section Error');
grid on;
legend('Uniform Mesh Error','Adaptive Mesh Error','Location','best');


% -----------------------
% 5. INTEGRATED ERROR ANALYSIS (Away from Singularity)
% -----------------------
% Build masks in normalized coordinates (consistent with mesh nodes)
R0_norm = R0_phys / phys_scale;   % proton radius in normalized units (same units as nodes)
mask_interior = @(nodes) (nodes(:,1) > 0.1*L) & (nodes(:,1) < 0.9*L) & ...
                         (nodes(:,2) > 0.1*L) & (nodes(:,2) < 0.9*L);

mask_radius   = @(nodes) sqrt((nodes(:,1)-center(1)).^2 + (nodes(:,2)-center(2)).^2) > 5*R0_norm;

mask_uni  = mask_interior(nodes_uni)  & mask_radius(nodes_uni);
mask_amr  = mask_interior(nodes_final)& mask_radius(nodes_final);

% compute errors (RMS and Linf) on masked nodes
err_nodal_uni = phi_uni(mask_uni) - phi_exact_uni(mask_uni);
err_nodal_amr = phi_plot_amr(mask_amr) - phi_exact_final(mask_amr);

L2_err_uni  = sqrt(mean(err_nodal_uni.^2));
Linf_err_uni= max(abs(err_nodal_uni));

L2_err_amr  = sqrt(mean(err_nodal_amr.^2));
Linf_err_amr= max(abs(err_nodal_amr));

fprintf('\n--- 5. Integrated Error Analysis (away from singularity) ---\n');
fprintf('Uniform mesh nodes = %d:   L2 (RMS) = %.3e   Linf = %.3e   (masked nodes = %d)\n', ...
        size(nodes_uni,1), L2_err_uni, Linf_err_uni, nnz(mask_uni));
fprintf('Adaptive mesh nodes = %d:  L2 (RMS) = %.3e   Linf = %.3e   (masked nodes = %d)\n', ...
        size(nodes_final,1), L2_err_amr, Linf_err_amr, nnz(mask_amr));

% relative numbers
if L2_err_uni > 0
    rel_red = (L2_err_uni - L2_err_amr) / L2_err_uni;
    fprintf('Adaptive vs Uniform: absolute L2 reduction = %.3e, relative = %.2f%%\n', L2_err_uni - L2_err_amr, rel_red*100);
    fprintf('Error ratio (L2 adaptive / L2 uniform) = %.3f\n', L2_err_amr / L2_err_uni);
end


%% =====================================================================
% HELPER FUNCTIONS (local)
% Note: local functions are defined below the script body
% =====================================================================

function [nodes, elem] = generate_initial_mesh(L, Nx, Ny)
    x_coords = linspace(0, L, Nx+1);
    y_coords = linspace(0, L, Ny+1);
    [Xg, Yg] = meshgrid(x_coords, y_coords);
    nodes = [Xg(:), Yg(:)];
    elem = zeros(2*Nx*Ny, 3);
    N_x = Nx + 1;
    idx = 1;
    for j = 1:Ny
        for i = 1:Nx
            n1 = (j-1)*N_x + i;
            n2 = n1 + 1;
            n3 = n1 + N_x;
            n4 = n2 + N_x;
            elem(idx, :) = [n1 n2 n3]; idx = idx + 1;
            elem(idx, :) = [n2 n4 n3]; idx = idx + 1;
        end
    end
    elem = elem(1:idx-1, :);
end

function [A, b, dN_dx_list, dN_dy_list, Area_list] = assemble_system(nodes, elem, N, Ne, center, rho_func)
    A = sparse(N, N);
    b = zeros(N, 1);
    dN_dx_list = zeros(3, Ne);
    dN_dy_list = zeros(3, Ne);
    Area_list  = zeros(Ne, 1);

    for e = 1:Ne
        ids = elem(e, :);
        x = nodes(ids, 1);
        y = nodes(ids, 2);

        Area = 0.5 * abs( x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2)) );
        Area_list(e) = Area;

        if Area < eps
            continue;
        end

        dN_dx = [y(2)-y(3); y(3)-y(1); y(1)-y(2)] / (2*Area);
        dN_dy = [x(3)-x(2); x(1)-x(3); x(2)-x(1)] / (2*Area);
        dN_dx_list(:, e) = dN_dx;
        dN_dy_list(:, e) = dN_dy;

        Ke = Area * (dN_dx * dN_dx' + dN_dy * dN_dy');
        A(ids, ids) = A(ids, ids) + Ke;

        xc = mean(x);
        yc = mean(y);
        rho = rho_func(xc, yc);
        be = (Area * rho / 3) * ones(3,1);
        b(ids) = b(ids) + be;
    end
end

function eta = error_estimation(nodes, elem, phi, dN_dx_list, dN_dy_list, Area_list, N, Ne)
    % Edge-jump estimator (residual-like). Returns per-element indicator.
    phi_x = zeros(Ne, 1);
    phi_y = zeros(Ne, 1);

    for e = 1:Ne
        phi_e = phi(elem(e, :));
        phi_x(e) = dN_dx_list(:, e)' * phi_e;
        phi_y(e) = dN_dy_list(:, e)' * phi_e;
    end

    TR = triangulation(elem, nodes);
    E_adj = edges(TR);
    T_adj = edgeAttachments(TR, E_adj);

    eta_sq = zeros(Ne, 1);

    for k = 1:size(E_adj, 1)
        T_ids = T_adj{k};
        if length(T_ids) ~= 2
            continue;
        end

        e1 = T_ids(1);
        e2 = T_ids(2);

        grad_jump_sq = (phi_x(e1) - phi_x(e2))^2 + (phi_y(e1) - phi_y(e2))^2;
        jump_mag = sqrt(grad_jump_sq);

        n_ids = E_adj(k, :);
        h_E = norm(nodes(n_ids(1), :) - nodes(n_ids(2), :));

        indicator = jump_mag * h_E;

        eta_sq(e1) = eta_sq(e1) + indicator^2;
        eta_sq(e2) = eta_sq(e2) + indicator^2;
    end

    eta = sqrt(eta_sq);
    if max(eta) > 0
        eta = eta / max(eta);
    end
end

function marked_elements = mark_elements(eta, theta)
    eta_sq = eta.^2;
    total_error_sq = sum(eta_sq);
    if total_error_sq == 0
        marked_elements = [];
        return;
    end
    target_error_sq = theta * total_error_sq;

    [~, sorted_indices] = sort(eta_sq, 'descend');

    cumulative_error_sq = 0;
    marked_count = 0;

    for k = 1:length(sorted_indices)
        cumulative_error_sq = cumulative_error_sq + eta_sq(sorted_indices(k));
        marked_count = marked_count + 1;
        if cumulative_error_sq >= target_error_sq
            break;
        end
    end

    marked_elements = sorted_indices(1:marked_count);
end

function marked_closed = closure_neighbours(marked_elements, elem)
    Ne = size(elem,1);
    neighbors = build_neighbors(elem);  % 1-ring neighbors
    marked = false(Ne,1);
    if ~isempty(marked_elements)
        marked(marked_elements) = true;
    end

    % Expand to 1-ring only
    current = find(marked);
    for e = current'
        marked(neighbors{e}) = true;
    end

    marked_closed = find(marked);
end

function [new_nodes, new_elem, edge_mid_map] = refine_mesh(nodes, elem, marked_elements, edge_mid_map)
    if nargin < 4 || isempty(edge_mid_map)
        edge_mid_map = containers.Map('KeyType','char','ValueType','int32');
    end

    N_old   = size(nodes, 1);
    Ne_old  = size(elem, 1);

    new_nodes = nodes;
    next_node_id = N_old + 1;

    new_elem = zeros(4*length(marked_elements) + (Ne_old-length(marked_elements)), 3);
    idx = 1;

    marked_mask = false(Ne_old,1);
    if ~isempty(marked_elements)
        marked_mask(marked_elements) = true;
    end

    for e = 1:Ne_old
        tri = elem(e,:);
        n1 = tri(1); n2 = tri(2); n3 = tri(3);

        if ~marked_mask(e)
            new_elem(idx,:) = tri;
            idx = idx + 1;
            continue;
        end

        edges = [n1 n2; n2 n3; n3 n1];
        mids  = zeros(3,1);

        for k = 1:3
            a = edges(k,1);
            b = edges(k,2);
            key = sprintf('%d-%d', min(a,b), max(a,b));

            if isKey(edge_mid_map, key)
                mids(k) = edge_mid_map(key);
            else
                new_nodes(next_node_id,:) = 0.5*(nodes(a,:) + nodes(b,:));
                edge_mid_map(key) = next_node_id;
                mids(k) = next_node_id;
                next_node_id = next_node_id + 1;
            end
        end

        m12 = mids(1); m23 = mids(2); m31 = mids(3);

        new_elem(idx,:) = [n1, m12, m31]; idx = idx+1;
        new_elem(idx,:) = [n2, m23, m12]; idx = idx+1;
        new_elem(idx,:) = [n3, m31, m23]; idx = idx+1;
        new_elem(idx,:) = [m12, m23, m31]; idx = idx+1;
    end

    new_elem = new_elem(1:idx-1,:);

    maxElemIndex = max(new_elem(:));
    nNodes = size(new_nodes,1);
    if maxElemIndex > nNodes
        error('REFINE_MESH: invalid mesh: element references node %d but only %d nodes exist', maxElemIndex, nNodes);
    end
end

function neighbors = build_neighbors(elem)
    Ne = size(elem,1);
    edge_map = containers.Map('KeyType','char','ValueType','any');

    for e = 1:Ne
        tri = elem(e,:);
        edges = [tri([1 2]); tri([2 3]); tri([3 1])];
        for k = 1:3
            edge = sort(edges(k,:));
            key = sprintf('%d-%d', edge(1), edge(2));
            if isKey(edge_map,key)
                edge_map(key) = [edge_map(key), e];
            else
                edge_map(key) = e;
            end
        end
    end

    neighbors = cell(Ne,1);
    for e = 1:Ne
        neighbors{e} = [];
    end

    vals = values(edge_map);
    for k = 1:length(vals)
        elems_sharing_edge = vals{k};
        if numel(elems_sharing_edge) < 2
            continue;
        end
        for i = 1:numel(elems_sharing_edge)
            e1 = elems_sharing_edge(i);
            others = setdiff(elems_sharing_edge, e1);
            neighbors{e1} = unique([neighbors{e1}, others]);
        end
    end
end

function Vq = triinterp(TR, V, Pq)
    % TRIINTERP - barycentric interpolation on an unstructured triangulation TR.
    nNodesTR = size(TR.Points,1);
    if length(V) ~= nNodesTR
        error('TRIINTERP: length(V) = %d incompatible with TR (%d nodes). Use nodal vector matching TR.', length(V), nNodesTR);
    end

    elem_idx = pointLocation(TR, Pq);
    valid_q = ~isnan(elem_idx);

    Vq = nan(size(Pq, 1), 1);

    for i = find(valid_q)'
        e = elem_idx(i);
        ids = TR.ConnectivityList(e, :);

        if any(ids > nNodesTR) || any(ids < 1)
            error('TRIINTERP: triangulation refers to nonexisting node. MaxID=%d, Nnodes=%d', max(ids), nNodesTR);
        end

        nodes_e = TR.Points(ids, :);

        x = nodes_e(:, 1); y = nodes_e(:, 2);
        x_q = Pq(i, 1); y_q = Pq(i, 2);

        Area_e = 0.5 * abs(x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2)));
        if Area_e == 0
            Vq(i) = NaN;
            continue;
        end

        L1 = (x_q*(y(2)-y(3)) + x(2)*(y(3)-y_q) + x(3)*(y_q-y(2))) / (2*Area_e);
        L2 = (x_q*(y(3)-y(1)) + x(3)*(y(1)-y_q) + x(1)*(y_q-y(3))) / (2*Area_e);
        L3 = 1 - L1 - L2;

        Vq(i) = L1 * V(ids(1)) + L2 * V(ids(2)) + L3 * V(ids(3));
    end
end
